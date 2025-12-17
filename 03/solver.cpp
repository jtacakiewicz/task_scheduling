#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <numeric>

using namespace std;

const int NUM_MACHINES = 5;

struct JobData {
    int id;
    vector<int> processing_times; 
    int due_date;
};

pair<long long, vector<int>> calculate_total_tardiness(
    const vector<int>& sequence, 
    const vector<JobData>& all_jobs, 
    const vector<vector<int>>& setup_times) 
{
    int n = sequence.size();
    vector<vector<int>> C(NUM_MACHINES, vector<int>(n, 0)); 
    long long total_tardiness = 0;
    vector<int> completion_on_M5(n, 0);

    for (int j = 0; j < n; ++j) {
        int job_idx = sequence[j];
        const auto& job = all_jobs[job_idx];

        for (int m = 0; m < NUM_MACHINES; ++m) {
            int start_time;
            if (m == 0) {
                int prev_C_M1 = (j == 0) ? 0 : C[0][j-1]; 
                int prev_job_idx = (j == 0) ? -1 : sequence[j-1]; // -1 dla pierwszego zadania
                int setup_time = (j == 0) ? 0 : setup_times[prev_job_idx][job_idx];
                
                start_time = prev_C_M1 + setup_time;
            } 
            else { 
                start_time = C[m-1][j];
            }

            if (j > 0) {
                int prev_job_idx = sequence[j-1];
                int prev_C_Mm = C[m][j-1];
                int setup_time = setup_times[prev_job_idx][job_idx];

                start_time = max(start_time, prev_C_Mm + setup_time);
            }
            
            C[m][j] = start_time + job.processing_times[m];
        }

        int finish_time = C[NUM_MACHINES - 1][j];
        completion_on_M5[j] = finish_time;
        int tardiness = max(0, finish_time - job.due_date);
        total_tardiness += tardiness;
    }

    return {total_tardiness, completion_on_M5};
}


class FlowShopSolver {
private:
    int n;
    vector<JobData> all_jobs;
    vector<vector<int>> setup_times;
    const int MAX_ITERATIONS = 10000;
    const int TABU_TENURE = 10;

    struct Move {
        int i, j;
        long long cost;
        vector<int> new_sequence;
    };

    pair<long long, vector<int>> tabu_search() {
        vector<int> current_sequence(n);
        iota(current_sequence.begin(), current_sequence.end(), 0);
        long long current_cost = calculate_total_tardiness(current_sequence, all_jobs, setup_times).first;

        vector<int> best_sequence = current_sequence;
        long long best_cost = current_cost;

        vector<vector<int>> tabu_list(n, vector<int>(n, 0)); 
        int iteration = 0;

        while (iteration < MAX_ITERATIONS) {
            Move best_move_in_neighborhood = {-1, -1, -1, {}};
            
            for (int i = 0; i < n; ++i) {
                for (int j = i + 1; j < n; ++j) {
                    
                    vector<int> neighbor_sequence = current_sequence;
                    swap(neighbor_sequence[i], neighbor_sequence[j]);
                    
                    long long neighbor_cost = calculate_total_tardiness(neighbor_sequence, all_jobs, setup_times).first;
                    int job_i = current_sequence[i];
                    int job_j = current_sequence[j]; 
                    
                    bool is_tabu = (tabu_list[job_i][job_j] > iteration || tabu_list[job_j][job_i] > iteration);

                    if (neighbor_cost < best_cost || !is_tabu) {
                        
                        if (best_move_in_neighborhood.cost == -1 || neighbor_cost < best_move_in_neighborhood.cost) {
                            best_move_in_neighborhood = {i, j, neighbor_cost, neighbor_sequence};
                        }
                    }
                }
            }

            if (best_move_in_neighborhood.cost == -1) {
                break; 
            }

            current_sequence = best_move_in_neighborhood.new_sequence;
            current_cost = best_move_in_neighborhood.cost;

            int job_i_moved_to_j = current_sequence[best_move_in_neighborhood.j];
            int job_j_moved_to_i = current_sequence[best_move_in_neighborhood.i];
            
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
        if (!file.is_open()) {
            cerr << "Błąd: Nie można otworzyć pliku " << filename << endl;
            return false;
        }

        string line;
        if (getline(file, line)) {
            stringstream ss(line);
            ss >> n;
            if (n <= 0) {
                 cerr << "Błąd: Nieprawidłowa liczba zadań n." << endl;
                 return false;
            }
        } else {
            cerr << "Błąd: Pusty plik lub błąd odczytu n." << endl;
            return false;
        }

        all_jobs.resize(n);
        setup_times.resize(n, vector<int>(n));
        for (int i = 0; i < n; ++i) {
            if (getline(file, line)) {
                stringstream ss(line);
                all_jobs[i].id = i + 1; // 1-bazowy ID
                all_jobs[i].processing_times.resize(NUM_MACHINES);
                for (int m = 0; m < NUM_MACHINES; ++m) {
                    ss >> all_jobs[i].processing_times[m];
                }
                ss >> all_jobs[i].due_date;
            } else {
                cerr << "Błąd: Oczekiwano " << n << " linii z czasami przetwarzania i terminami." << endl;
                return false;
            }
        }

        for (int i = 0; i < n; ++i) {
            if (getline(file, line)) {
                stringstream ss(line);
                for (int j = 0; j < n; ++j) {
                    ss >> setup_times[i][j];
                }
            } else {
                cerr << "Błąd: Oczekiwano " << n << " linii dla macierzy Sij." << endl;
                return false;
            }
        }
        
        file.close();
        return true;
    }

    void solve_and_save(const string& output_filename) {
        if (n == 0) {
            cerr << "Błąd: Brak danych do rozwiązania." << endl;
            return;
        }

        pair<long long, vector<int>> result = tabu_search();
        long long best_cost = result.first;
        const vector<int>& best_sequence_indices = result.second;

        ofstream outfile(output_filename);
        if (!outfile.is_open()) {
            cerr << "Błąd: Nie można otworzyć pliku do zapisu " << output_filename << endl;
            return;
        }

        outfile << best_cost << endl;

        for (size_t i = 0; i < best_sequence_indices.size(); ++i) {
            outfile << all_jobs[best_sequence_indices[i]].id << (i == best_sequence_indices.size() - 1 ? "" : " ");
        }
        outfile << endl;

        outfile.close();
        cout << "Rozwiązanie zapisane do pliku: " << output_filename << endl;
    }
};


int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Użycie: " << argv[0] << " <plik_wejsciowy> <plik_wynikowy>" << endl;
        return 1;
    }

    string input_file = argv[1];
    string output_file = argv[2];

    FlowShopSolver solver;

    if (solver.load_data(input_file)) {
        solver.solve_and_save(output_file);
    }

    return 0;
}
