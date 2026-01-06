#include <algorithm>
#include <chrono> // Dodano dla obsługi czasu
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
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
    pair<unsigned long long, vector<unsigned long long>> calculate_total_tardiness(const vector<int> &sequence) {
        int seq_size = sequence.size();
        vector<vector<unsigned long long>> C(NUM_MACHINES, vector<unsigned long long>(seq_size, 0));
        unsigned long long total_tardiness = 0;

        for (int j = 0; j < seq_size; ++j) {
            int current_job_idx = sequence[j];
            int prev_job_idx = (j == 0) ? -1 : sequence[j - 1];

            for (int m = 0; m < NUM_MACHINES; ++m) {
                unsigned long long machine_ready = 0;
                if (j > 0) {
                    machine_ready = C[m][j - 1] + setup_times[prev_job_idx][current_job_idx];
                }
                unsigned long long prev_stage_ready = (m == 0) ? 0 : C[m - 1][j];

                C[m][j] = max(machine_ready, prev_stage_ready) + all_jobs[current_job_idx].processing_times[m];
            }
            
            unsigned long long finish_time = C[NUM_MACHINES - 1][j];
            unsigned long long due = all_jobs[current_job_idx].due_date;
            if (finish_time > due) total_tardiness += (finish_time - due);
        }
        return {total_tardiness, C[NUM_MACHINES - 1]};
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
        unsigned long long current_D = calculate_total_tardiness(seq).first;
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
                    unsigned long long new_D = calculate_total_tardiness(seq).first;

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

        for (int i = 0; i < 100; i++) {
            double a = ((float)i) / 100;
            auto seq = get_setup_weighted_sequence(a);
            if(i == 0) {
                seq=get_edd_sequence();
            }
            unsigned long long d = calculate_total_tardiness(seq).first;
            if (best_overall_seq.empty() || d < best_overall_D) {
                best_overall_D = d;
                best_overall_seq = seq;
            }
        }

        auto result = local_search(best_overall_seq, (double)time_limit_sec * 0.95, start_time);
        
        return result;
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
