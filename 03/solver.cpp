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
    pair<unsigned long long, vector<unsigned long long>> calculate_total_tardiness(
        const vector<int> &sequence) {
        int n = sequence.size();
        vector<vector<unsigned long long>> C(NUM_MACHINES,
                                             vector<unsigned long long>(n, 0));
        unsigned long long total_tardiness = 0;
        vector<unsigned long long> completion_on_M5(n, 0);

        for (int j = 0; j < n; ++j) {
            int job_idx = sequence[j];
            const auto &job = all_jobs[job_idx];

            for (int m = 0; m < NUM_MACHINES; ++m) {
                unsigned long long start_time;
                if (m == 0) {
                    unsigned long long prev_C_M1 = (j == 0) ? 0 : C[0][j - 1];
                    int prev_job_idx = (j == 0) ? -1 : sequence[j - 1];
                    unsigned long long setup_time =
                        (j == 0) ? 0 : setup_times[prev_job_idx][job_idx];
                    start_time = prev_C_M1 + setup_time;
                } else {
                    start_time = C[m - 1][j];
                }

                if (j > 0) {
                    int prev_job_idx = sequence[j - 1];
                    unsigned long long prev_C_Mm = C[m][j - 1];
                    unsigned long long setup_time = setup_times[prev_job_idx][job_idx];
                    start_time = max(start_time, prev_C_Mm + setup_time);
                }

                C[m][j] = start_time + job.processing_times[m];
            }

            unsigned long long finish_time = C[NUM_MACHINES - 1][j];
            completion_on_M5[j] = finish_time;
            unsigned long long tardiness =
                (finish_time > job.due_date) ? (finish_time - job.due_date) : 0;
            total_tardiness += tardiness;
        }
        return {total_tardiness, completion_on_M5};
    }

    pair<unsigned long long, vector<int>> select_based_on_setup(double a) {
        vector<int> unscheduled(n);
        iota(unscheduled.begin(), unscheduled.end(), 0);

        vector<int> scheduled;

        while (!unscheduled.empty()) {
            int last = scheduled.empty() ? -1 : scheduled.back();
            int best_job = -1;
            double best_val = 1e18;

            for (int j : unscheduled) {
                unsigned long long setup =
                    (last == -1 ? 0 : setup_times[last][j]);
                double value = all_jobs[j].due_date + a * setup;

                if (value < best_val) {
                    best_val = value;
                    best_job = j;
                }
            }

            scheduled.push_back(best_job);
            unscheduled.erase(find(unscheduled.begin(),
                                   unscheduled.end(), best_job));
        }

        unsigned long long D = calculate_total_tardiness(scheduled).first;
        return make_pair(D, scheduled);
    }
    pair<unsigned long long, vector<int>> greed(vector<int> seq,
                         unsigned long long best_D,
                         double time_limit_sec) {
        auto start = chrono::steady_clock::now();
        unsigned long long current_D = best_D;
        bool improved = true;

        while (improved) {
            improved = false;

            for (int i = 0; i < n; ++i) {
                for (int j = i + 1; j < n; ++j) {
                    auto now = chrono::steady_clock::now();
                    double elapsed =
                        chrono::duration<double>(now - start).count();
                    if (elapsed > time_limit_sec)
                        return make_pair(current_D, seq);

                    vector<int> new_seq = seq;
                    swap(new_seq[i], new_seq[j]);

                    unsigned long long D = calculate_total_tardiness(new_seq).first;
                    if (D < current_D) {
                        seq = new_seq;
                        current_D = D;
                        improved = true;
                        goto NEXT_ITER;
                    }
                }
            }
        NEXT_ITER:;
        }
        return make_pair(current_D, seq);
    };

    pair<unsigned long long, vector<int>> solve(int time_limit_sec) {

        auto try_default_but_different_setup = [&]() {
            vector<int> seq(n);
            iota(seq.begin(), seq.end(), 0);
            if (n >= 2) swap(seq[0], seq[1]);
            return make_pair(seq, calculate_total_tardiness(seq).first);
        };
        auto start_time = chrono::steady_clock::now();

        vector<double> steps;
        for (int i = 0; i < 10; ++i)
            steps.push_back(i / 10.0);

        auto best = select_based_on_setup(steps[0]);
        vector<int> best_jobs = best.second;
        unsigned long long best_D = best.first;

        for (size_t i = 1; i < steps.size(); ++i) {
            auto res = select_based_on_setup(steps[i]);
            if (res.first < best_D) {
                best_D = res.first;
                best_jobs = res.second;
            }
        }

        auto def = try_default_but_different_setup();
        if (def.second < best_D) {
            best_D = def.second;
            best_jobs = def.first;
        }

        auto potential = select_based_on_setup(0.1);
        if (potential.first < best_D) {
            best_D = potential.first;
            best_jobs = potential.second;
        }

        auto greedy =
            greed(best_jobs, best_D, 0.9 * time_limit_sec);

        if (greedy.first < best_D) {
            best_D = greedy.first;
            best_jobs = greedy.second;
        }

        return { best_D, best_jobs };
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
