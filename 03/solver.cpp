#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <numeric>
#include <chrono> // Dodano dla obsługi czasu

using namespace std;

const int NUM_MACHINES = 5;

struct JobData {
    int id;
    vector<unsigned long long> processing_times; 
    unsigned long long due_date;
};

pair<unsigned long long, vector<unsigned long long>> calculate_total_tardiness(
    const vector<int>& sequence, 
    const vector<JobData>& all_jobs, 
    const vector<vector<unsigned long long>>& setup_times) 
{
    int n = sequence.size();
    vector<vector<unsigned long long>> C(NUM_MACHINES, vector<unsigned long long>(n, 0)); 
    unsigned long long total_tardiness = 0;
    vector<unsigned long long> completion_on_M5(n, 0);

    for (int j = 0; j < n; ++j) {
        int job_idx = sequence[j];
        const auto& job = all_jobs[job_idx];

        for (int m = 0; m < NUM_MACHINES; ++m) {
            unsigned long long start_time;
            if (m == 0) {
                unsigned long long prev_C_M1 = (j == 0) ? 0 : C[0][j-1]; 
                int prev_job_idx = (j == 0) ? -1 : sequence[j-1];
                unsigned long long setup_time = (j == 0) ? 0 : setup_times[prev_job_idx][job_idx];
                start_time = prev_C_M1 + setup_time;
            } 
            else { 
                start_time = C[m-1][j];
            }

            if (j > 0) {
                int prev_job_idx = sequence[j-1];
                unsigned long long prev_C_Mm = C[m][j-1];
                unsigned long long setup_time = setup_times[prev_job_idx][job_idx];
                start_time = max(start_time, prev_C_Mm + setup_time);
            }
            
            C[m][j] = start_time + job.processing_times[m];
        }

        unsigned long long finish_time = C[NUM_MACHINES - 1][j];
        completion_on_M5[j] = finish_time;
        unsigned long long tardiness = (finish_time > job.due_date) ? (finish_time - job.due_date) : 0;
        total_tardiness += tardiness;
    }
    return {total_tardiness, completion_on_M5};
}

class FlowShopSolver {
private:
    int n;
    vector<JobData> all_jobs;
    vector<vector<unsigned long long>> setup_times;
    const int MAX_ITERATIONS = 1000000; // Zwiększono limit iteracji, bo teraz ogranicza nas czas
    const int TABU_TENURE = 10;

    struct Move {
        int i, j;
        unsigned long long cost;
        vector<int> new_sequence;
    };

    // Dodano parametr time_limit_sec
    pair<unsigned long long, vector<int>> tabu_search(int time_limit_sec) {
        auto start_wall_time = chrono::steady_clock::now();
        
        vector<int> current_sequence(n);
        iota(current_sequence.begin(), current_sequence.end(), 0);
        unsigned long long current_cost = calculate_total_tardiness(current_sequence, all_jobs, setup_times).first;

        vector<int> best_sequence = current_sequence;
        unsigned long long best_cost = current_cost;

        vector<vector<int>> tabu_list(n, vector<int>(n, 0)); 
        int iteration = 0;
        const unsigned long long ULL_MAX_VAL = ~0ULL;

        while (iteration < MAX_ITERATIONS) {
            // SPRAWDZENIE LIMITU CZASU
            auto current_now = chrono::steady_clock::now();
            auto elapsed = chrono::duration_cast<chrono::milliseconds>(current_now - start_wall_time).count();
            if (elapsed+500 >= time_limit_sec*1000) {
                break;
            }

            Move best_move_in_neighborhood = {-1, -1, ULL_MAX_VAL, {}};
            
            for (int i = 0; i < n; ++i) {
                for (int j = i + 1; j < n; ++j) {
                    vector<int> neighbor_sequence = current_sequence;
                    swap(neighbor_sequence[i], neighbor_sequence[j]);
                    
                    unsigned long long neighbor_cost = calculate_total_tardiness(neighbor_sequence, all_jobs, setup_times).first;
                    int job_i = current_sequence[i];
                    int job_j = current_sequence[j]; 
                    
                    bool is_tabu = (tabu_list[job_i][job_j] > iteration || tabu_list[job_j][job_i] > iteration);

                    if (neighbor_cost < best_cost || !is_tabu) {
                        if (neighbor_cost < best_move_in_neighborhood.cost) {
                            best_move_in_neighborhood = {i, j, neighbor_cost, neighbor_sequence};
                        }
                    }
                }
            }

            if (best_move_in_neighborhood.cost == ULL_MAX_VAL) break;

            current_sequence = best_move_in_neighborhood.new_sequence;
            current_cost = best_move_in_neighborhood.cost;

            int job1 = current_sequence[best_move_in_neighborhood.i];
            int job2 = current_sequence[best_move_in_neighborhood.j];
            tabu_list[job1][job2] = iteration + TABU_TENURE;
            tabu_list[job2][job1] = iteration + TABU_TENURE;

            if (current_cost < best_cost) {
                best_cost = current_cost;
                best_sequence = current_sequence;
            }
            iteration++;
        }
        return {best_cost, best_sequence};
    }

public:
    FlowShopSolver() : n(0) {}

    bool load_data(const string& filename) {
        ifstream file(filename);
        if (!file.is_open()) return false;
        string line;
        if (getline(file, line)) {
            stringstream ss(line);
            ss >> n;
        } else return false;

        all_jobs.resize(n);
        setup_times.assign(n, vector<unsigned long long>(n));
        for (int i = 0; i < n; ++i) {
            if (getline(file, line)) {
                stringstream ss(line);
                all_jobs[i].id = i + 1;
                all_jobs[i].processing_times.resize(NUM_MACHINES);
                for (int m = 0; m < NUM_MACHINES; ++m) ss >> all_jobs[i].processing_times[m];
                ss >> all_jobs[i].due_date;
            }
        }
        for (int i = 0; i < n; ++i) {
            if (getline(file, line)) {
                stringstream ss(line);
                for (int j = 0; j < n; ++j) ss >> setup_times[i][j];
            }
        }
        return true;
    }

    void solve_and_save(const string& output_filename, int time_limit_sec) {
        if (n == 0) return;

        auto result = tabu_search(time_limit_sec);
        unsigned long long best_cost = result.first;
        const vector<int>& best_seq = result.second;

        ofstream outfile(output_filename);
        if (outfile.is_open()) {
            outfile << best_cost << endl;
            for (size_t i = 0; i < best_seq.size(); ++i) {
                outfile << all_jobs[best_seq[i]].id << (i == best_seq.size() - 1 ? "" : " ");
            }
            outfile << endl;
        }
    }
};

int main(int argc, char* argv[]) {
    if (argc != 4) {
        cerr << "Użycie: " << argv[0] << " <wejście> <wyjście> <limit_sekundy>" << endl;
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
