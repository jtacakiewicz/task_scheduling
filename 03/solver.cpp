#include <algorithm>
#include <chrono> // Dodano dla obsługi czasu
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <vector>

using namespace std;

const int NUM_MACHINES = 5;

struct JobData {
    int id;
    vector<unsigned long long> processing_times;
    unsigned long long due_date;
};

class FlowShopSolver {
private:
    int n;
    vector<JobData> all_jobs;
    vector<vector<unsigned long long>> setup_times;
    const int MAX_ITERATIONS = 1000000;
    mt19937 rng{2137};

    unsigned long long machine_finish_times[NUM_MACHINES];

    pair<unsigned long long, vector<unsigned long long>> get_tardiness(const vector<int> &sequence) {
        int seq_size = sequence.size();
        vector<unsigned long long> last_machine_C(seq_size);
        unsigned long long total_tardiness = 0;

        for (int m = 0; m < NUM_MACHINES; ++m) machine_finish_times[m] = 0;

        int prev_job_idx = -1;

        for (int j = 0; j < seq_size; ++j) {
            int current_job_idx = sequence[j];

            unsigned long long setup = (prev_job_idx == -1) ? 0 : setup_times[prev_job_idx][current_job_idx];
            machine_finish_times[0] = machine_finish_times[0] + setup + all_jobs[current_job_idx].processing_times[0];

            for (int m = 1; m < NUM_MACHINES; ++m) {
                unsigned long long ready_from_prev_job = machine_finish_times[m] + setup;
                unsigned long long ready_from_prev_stage = machine_finish_times[m-1];

                machine_finish_times[m] = (ready_from_prev_job > ready_from_prev_stage ? 
                    ready_from_prev_job : ready_from_prev_stage) 
                    + all_jobs[current_job_idx].processing_times[m];
            }

            unsigned long long finish_time = machine_finish_times[NUM_MACHINES - 1];
            last_machine_C[j] = finish_time;

            unsigned long long due = all_jobs[current_job_idx].due_date;
            if (finish_time > due) {
                total_tardiness += (finish_time - due);
            }

            prev_job_idx = current_job_idx;
        }

        return {total_tardiness, last_machine_C};
    }

    void insert_best(vector<int>& seq, int job) {
        int best_pos = 0;
        unsigned long long min_d = -1; // max

        for (int i = 0; i <= seq.size(); ++i) {
            seq.insert(seq.begin() + i, job);
            unsigned long long d = get_tardiness(seq).first;
            if (d < min_d) {
                min_d = d;
                best_pos = i;
            }
            seq.erase(seq.begin() + i);
        }
        seq.insert(seq.begin() + best_pos, job);
    }
    void adjacent_local_search(vector<int>& seq, unsigned long long& current_d) {
        bool improved = true;
        while (improved) {
            improved = false;
            for (int i = 0; i < (int)seq.size() - 1; ++i) {
                swap(seq[i], seq[i+1]);
                unsigned long long new_d = get_tardiness(seq).first;
                
                if (new_d < current_d) {
                    current_d = new_d;
                    improved = true;
                } else {
                    swap(seq[i], seq[i+1]);
                }
            }
        }
    }
    pair<unsigned long long, vector<int>> solve_small(int time_limit_sec, vector<int> current_seq) {
        auto start_time = chrono::steady_clock::now();
        double time_limit = (double)time_limit_sec;

        unsigned long long current_d = get_tardiness(current_seq).first;
        adjacent_local_search(current_seq, current_d);

        vector<int> best_seq = current_seq;
        unsigned long long best_d = current_d;

        // --- INICJALIZACJA TEMPERATURY ---
        // Dobra praktyka: T0 to ok. 20-50% średniego kosztu (opóźnienia) 
        // Możesz też ustawić stałą wartość, np. 100.0, i dostosować ją testowo.
        double T0 = max(1.0, 0.3 * current_d / current_seq.size()); 
        double temperature = T0;

        while (true) {
            auto now = chrono::steady_clock::now();
            double elapsed = chrono::duration<double>(now - start_time).count();
            if (elapsed >= time_limit) break;

            temperature = T0 * (1.0 - (elapsed / time_limit));
            if (temperature < 0.01) temperature = 0.01;

            vector<int> temp_seq = current_seq;

            int d_size = 3 + (rng() % 3);
            vector<int> removed;
            for (int i = 0; i < d_size; ++i) {
                if (temp_seq.empty()) break;
                int idx = uniform_int_distribution<int>(0, (int)temp_seq.size() - 1)(rng);
                removed.push_back(temp_seq[idx]);
                temp_seq.erase(temp_seq.begin() + idx);
            }

            for (int job : removed) {
                insert_best(temp_seq, job);
            }

            unsigned long long temp_d = get_tardiness(temp_seq).first;

            if (temp_d < current_d) {
                current_seq = temp_seq;
                current_d = temp_d;
                if (current_d < best_d) {
                    best_d = current_d;
                    best_seq = current_seq;
                    adjacent_local_search(best_seq, best_d);
                }
            } else {
                double delta = (double)(temp_d - current_d);
                double prob = exp(-delta / temperature);

                if (uniform_real_distribution<double>(0, 1)(rng) < prob) {
                    current_seq = temp_seq;
                    current_d = temp_d;
                }
            }
        }

        return {best_d, best_seq};
    }
    struct Individual {
        vector<int> sequence;
        unsigned long long tardiness;

        bool operator<(const Individual& other) const {
            return tardiness < other.tardiness;
        }
    };
    int pop_size = 40;
    double mutation_rate = 0.15;
    double adj_search_rate = 0.05;

    // Operator OX (Order Crossover) - zachowuje relacje kolejności
    vector<int> crossover_ox(const vector<int>& p1, const vector<int>& p2) {
        int n = p1.size();
        vector<int> child(n, -1);

        // Wybieramy losowy przedział z pierwszego rodzica
        int start = rng() % n;
        int end = rng() % n;
        if (start > end) swap(start, end);

        for (int i = start; i <= end; ++i) {
            child[i] = p1[i];
        }

        // Wypełniamy resztę genami z p2 w kolejności, w jakiej w nim występują
        int p2_idx = 0;
        for (int i = 0; i < n; ++i) {
            if (i >= start && i <= end) continue;

            while (true) {
                int gene = p2[p2_idx++];
                // Sprawdź czy gen już jest w dziecku
                bool exists = false;
                for (int j = start; j <= end; ++j) {
                    if (child[j] == gene) { exists = true; break; }
                }
                if (!exists) {
                    child[i] = gene;
                    break;
                }
            }
        }
        return child;
    }
    pair<unsigned long long, vector<int>> solve_genetic(double time_limit, vector<int> start_seq) {
        auto start_time = chrono::steady_clock::now();
        vector<Individual> population;

        // 1. Inicjalizacja populacji (start_seq + losowe mutacje + całkowicie losowe)
        population.push_back({start_seq, get_tardiness(start_seq).first});
        while (population.size() < pop_size) {
            vector<int> seq = start_seq;
            if (rng() % 2 == 0) mutate(seq);
            else shuffle(seq.begin(), seq.end(), rng);
            population.push_back({seq, get_tardiness(seq).first});
        }

        while (chrono::duration<double>(chrono::steady_clock::now() - start_time).count() < time_limit) {
            sort(population.begin(), population.end());

            vector<Individual> next_gen;
            // Elitaryzm (zachowaj 2 najlepszych)
            next_gen.push_back(population[0]);
            next_gen.push_back(population[1]);

            while (next_gen.size() < pop_size) {
                // Selekcja turniejowa (wybierz 2 z 5 losowych)
                auto tournament = [&]() {
                    int best_idx = rng() % pop_size;
                    for (int i = 0; i < 3; ++i) {
                        int idx = rng() % pop_size;
                        if (population[idx].tardiness < population[best_idx].tardiness) best_idx = idx;
                    }
                    return population[best_idx].sequence;
                };

                vector<int> p1 = tournament();
                vector<int> p2 = tournament();

                vector<int> child_seq = crossover_ox(p1, p2);

                if (uniform_real_distribution<double>(0, 1)(rng) < mutation_rate) {
                    mutate(child_seq);
                }

                auto tardiness = get_tardiness(child_seq).first;
                if (uniform_real_distribution<double>(0, 1)(rng) < adj_search_rate) {
                    adjacent_local_search(child_seq, tardiness);
                }

                next_gen.push_back({child_seq, tardiness});
            }
            population = next_gen;
        }

        sort(population.begin(), population.end());
        return {population[0].tardiness, population[0].sequence};
    }

    void mutate(vector<int>& seq) {
        // Swap mutation: zamień dwa losowe zadania
        int i = rng() % seq.size();
        int j = rng() % seq.size();
        swap(seq[i], seq[j]);
    }

    vector<int> get_neh_sequence() {
        vector<pair<unsigned long long, int>> total_processing;
        for (int j = 0; j < n; ++j) {
            unsigned long long sum_p = 0;
            for (int m = 0; m < NUM_MACHINES; ++m) {
                sum_p += all_jobs[j].processing_times[m];
            }
            total_processing.push_back({sum_p, j});
        }

        // Sortujemy malejąco po sumarycznym czasie przetwarzania
        sort(total_processing.rbegin(), total_processing.rend());

        vector<int> seq;
        // 2. Iteracyjne wstawianie zadań (procedura NEH)
        for (int i = 0; i < n; ++i) {
            int job_to_insert = total_processing[i].second;

            int best_pos = 0;
            unsigned long long min_tardiness = -1; // Max value

            // Szukamy najlepszej pozycji w aktualnej sekwencji
            for (int pos = 0; pos <= (int)seq.size(); ++pos) {
                seq.insert(seq.begin() + pos, job_to_insert);
                unsigned long long current_tardiness = get_tardiness(seq).first;

                if (current_tardiness < min_tardiness) {
                    min_tardiness = current_tardiness;
                    best_pos = pos;
                }
                seq.erase(seq.begin() + pos);
            }

            // Wstawiamy zadanie na stałe w najlepsze miejsce
            seq.insert(seq.begin() + best_pos, job_to_insert);
        }

        return seq;
    }

    vector<int> get_edd_sequence() {
        vector<int> seq(n);
        iota(seq.begin(), seq.end(), 0);
        sort(seq.begin(), seq.end(), [&](int a, int b) {
            return all_jobs[a].due_date < all_jobs[b].due_date;
        });
        return seq;
    }

    vector<int> get_setup_weighted_sequence(double a) {
        vector<int> unscheduled(n);
        iota(unscheduled.begin(), unscheduled.end(), 0);
        vector<int> scheduled;

        while (!unscheduled.empty()) {
            int last = scheduled.empty() ? -1 : scheduled.back();
            int best_idx_in_unscheduled = 0;
            double best_val = 1e18;

            for (int i = 0; i < unscheduled.size(); ++i) {
                int j = unscheduled[i];
                unsigned long long setup = (last == -1 ? 0 : setup_times[last][j]);
                double value = all_jobs[j].due_date + a * setup;
                if (value < best_val) {
                    best_val = value;
                    best_idx_in_unscheduled = i;
                }
            }
            scheduled.push_back(unscheduled[best_idx_in_unscheduled]);
            unscheduled.erase(unscheduled.begin() + best_idx_in_unscheduled);
        }
        return scheduled;
    }

    pair<unsigned long long, vector<int>> local_search(vector<int> seq, double time_limit_sec, chrono::steady_clock::time_point start_time) {
        unsigned long long current_D = get_tardiness(seq).first;
        bool improved = true;

        while (improved) {
            improved = false;
            for (int i = 0; i < n - 1; ++i) {
                for (int j = i + 1; j < n; ++j) {
                    auto now = chrono::steady_clock::now();
                    if (chrono::duration<double>(now - start_time).count() >= time_limit_sec) {
                        return {current_D, seq};
                    }

                    swap(seq[i], seq[j]);
                    unsigned long long new_D = get_tardiness(seq).first;

                    if (new_D < current_D) {
                        current_D = new_D;
                        improved = true;
                        goto next_iteration; 
                    } else {
                        swap(seq[i], seq[j]);
                    }
                }
            }
            next_iteration:;
        }
        return {current_D, seq};
    }

    pair<unsigned long long, vector<int>> solve(int time_limit_sec) {
        auto start_time = chrono::steady_clock::now();
        
        unsigned long long best_overall_D = -1; // max value
        vector<int> best_overall_seq;

        for (int i = 0; i <= 101; i++) {
            double a = ((float)i) / 100;
            vector<int> seq;
            if(i == 100) {
                seq = get_edd_sequence();
            }else if(i == 101) {
                seq = get_neh_sequence();
            }else {
                seq = get_setup_weighted_sequence(a);
            }
            unsigned long long d = get_tardiness(seq).first;
            if (best_overall_seq.empty() || d < best_overall_D) {
                best_overall_D = d;
                best_overall_seq = seq;
            }
        }
        auto now = chrono::steady_clock::now();
        double elapsed = chrono::duration<double>(now - start_time).count();
        double time_left = (static_cast<double>(time_limit_sec) - elapsed);
        auto result = solve_small(time_left * 0.5, best_overall_seq);
        return solve_with_ants(time_left * 0.5, result.second);
    }

public:
    FlowShopSolver() : n(0) {}

    bool load_data(const string &filename) {
        ifstream file(filename);
        if (!file.is_open())
            return false;
        string line;
        if (getline(file, line)) {
            stringstream ss(line);
            ss >> n;
        } else {
            return false;
        }
        all_jobs.resize(n);
        setup_times.assign(n, vector<unsigned long long>(n));
        for (int i = 0; i < n; ++i) {
            if (getline(file, line)) {
                stringstream ss(line);
                all_jobs[i].id = i + 1;
                all_jobs[i].processing_times.resize(NUM_MACHINES);
                for (int m = 0; m < NUM_MACHINES; ++m)
                    ss >> all_jobs[i].processing_times[m];
                ss >> all_jobs[i].due_date;
            }
        }
        for (int i = 0; i < n; ++i) {
            if (getline(file, line)) {
                stringstream ss(line);
                for (int j = 0; j < n; ++j)
                    ss >> setup_times[i][j];
            }
        }
        return true;
    }

    void solve_and_save(const string &output_filename, int time_limit_sec) {
        if (n == 0)
            return;

        auto result = solve(time_limit_sec);
        unsigned long long best_cost = result.first;
        const vector<int> &best_seq = result.second;

        ofstream outfile(output_filename);
        if (outfile.is_open()) {
            outfile << best_cost << endl;
            for (size_t i = 0; i < best_seq.size(); ++i) {
                outfile << all_jobs[best_seq[i]].id
                    << (i == best_seq.size() - 1 ? "" : " ");
            }
            outfile << endl;
        }
    }
};

int main(int argc, char *argv[]) {
    if (argc != 4) {
        cerr << "Użycie: " << argv[0] << " <wejście> <wyjście> <limit_sekundy>"
            << endl;
        return 1;
    }

    string input_file = argv[1];
    string output_file = argv[2];
    int time_limit = stoi(argv[3]);

    FlowShopSolver solver;
    if (solver.load_data(input_file)) {
        solver.solve_and_save(output_file, time_limit);
    }

    return 0;
}
