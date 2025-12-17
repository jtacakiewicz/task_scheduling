#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <random>
#include <ctime>

using namespace std;

const long long NUM_MACHINES = 5;

struct ProblemData {
    size_t n;
    vector<vector<long long>> processing_times;
    vector<long long> due_dates;
    vector<vector<long long>> setup_times;
};

void read_instance_file(const string& filename, ProblemData& data) {
    ifstream infile(filename);
    if (!infile.is_open()) {
        throw runtime_error("Nie można otworzyć pliku instancji: " + filename);
    }

    if (!(infile >> data.n)) {
        throw runtime_error("Błąd wczytywania liczby zamówień (n). Plik: " + filename);
    }

    data.processing_times.resize(data.n, vector<long long>(NUM_MACHINES));
    data.due_dates.resize(data.n);

    for (size_t i = 0; i < data.n; ++i) {
        for (size_t m = 0; m < NUM_MACHINES; ++m) {
            if (!(infile >> data.processing_times[i][m])) {
                throw runtime_error("Błąd wczytywania p_ij.");
            }
        }
        if (!(infile >> data.due_dates[i])) {
            throw runtime_error("Błąd wczytywania terminu d_j.");
        }
    }

    data.setup_times.resize(data.n, vector<long long>(data.n));
    for (size_t i = 0; i < data.n; ++i) {
        for (size_t j = 0; j < data.n; ++j) {
            if (!(infile >> data.setup_times[i][j])) {
                throw runtime_error("Błąd wczytywania czasu przezbrajania S_ij.");
            }
        }
    }
    infile.close();
}
void read_result_file(const string& filename, size_t n, long long& expected_total_tardiness, vector<long long>& schedule) {
    ifstream outfile(filename);
    if (!outfile.is_open()) {
        throw runtime_error("Nie można otworzyć pliku wynikowego (out): " + filename);
    }

    if (!(outfile >> expected_total_tardiness)) {
        throw runtime_error("Błąd wczytywania sumy opóźnień (SigmaDj) z pliku: " + filename);
    }

    schedule.clear();
    string job_id_str;
    while (outfile >> job_id_str) {
        try {
            long long job_index_1based = stoi(job_id_str);
            if (job_index_1based < 1 || job_index_1based > n) {
                throw runtime_error("Indeks zadania poza zakresem (1 do " + to_string(n) + "): J" + to_string(job_index_1based));
            }
            schedule.push_back(job_index_1based - 1); 
        } catch (const invalid_argument&) {
            throw runtime_error("Błąd parsowania numeru zadania z: " + job_id_str);
        }
    }

    if (schedule.size() != n) {
        throw runtime_error("Błąd: Wczytano " + to_string(schedule.size()) + " zadań, oczekiwano " + to_string(n) + ".");
    }
    
    vector<long long> check_sequence = schedule;
    sort(check_sequence.begin(), check_sequence.end());
    for(long long i = 0; i < n; ++i) {
        if(check_sequence[i] != i) {
            throw runtime_error("Błąd: Sekwencja musi zawierać każde zadanie J1 do J" + to_string(n) + " dokładnie raz. (J" + to_string(i) + ")");
        }
    }

    outfile.close();
}
long long evaluate_schedule(const ProblemData& data, const vector<long long>& schedule) {
    if (schedule.size() != data.n) { return -1; }

    vector<vector<long long>> C(NUM_MACHINES, vector<long long>(data.n, 0));

    for (long long k = 0; k < data.n; ++k) {
        long long current_job_idx = schedule[k]; 
        long long previous_job_idx = (k > 0) ? schedule[k - 1] : -1; 

        for (long long m = 0; m < NUM_MACHINES; ++m) {
            
            long long prev_finish_time_on_M = (k == 0) ? 0 : C[m][k - 1];
            
            long long setup_time = 0;
            if (previous_job_idx != -1) {
                setup_time = data.setup_times[previous_job_idx][current_job_idx];
            }
            long long machine_ready_time = prev_finish_time_on_M + setup_time;

            long long job_available_time = (m == 0) ? 0 : C[m - 1][k];
            
            long long start_time = max(machine_ready_time, job_available_time);
            long long processing_time = data.processing_times[current_job_idx][m];
            C[m][k] = start_time + processing_time;
        } 
    } 

    long long total_tardiness = 0;
    
    for (long long k = 0; k < data.n; ++k) {
        long long job_idx = schedule[k]; 
        long long Cj = C[NUM_MACHINES - 1][k]; 
        long long dj = data.due_dates[job_idx]; 
        
        long long Dj = max(0LL, Cj - dj); 
        total_tardiness += Dj;
    }

    return total_tardiness;
}
void generate_instance(const string& filename, size_t n, long long P_DEV_MAX = 15, long long P_MAX = 50, long long S_MAX = 10, double WSP_TARDY = 1.2, long long U_TARGET = 420) {
    if (n <= 0) {
        throw invalid_argument("Liczba zadań (n) musi być większa od 0.");
    }
    
    ofstream outfile(filename);
    ofstream outfile_sol(filename + ".opt");
    if (!outfile.is_open()) {
        throw runtime_error("Nie można utworzyć pliku wyjściowego: " + filename);
    }

    mt19937 generator(static_cast<unsigned long long>(time(0)));
    uniform_int_distribution<long long> p_dev_dist(-P_DEV_MAX/2, P_DEV_MAX/2);
    uniform_int_distribution<long long> p_dist(1, P_MAX);
    uniform_int_distribution<long long> s_dist(0, S_MAX);
    
    vector<vector<long long>> processing_times(n, vector<long long>(NUM_MACHINES));
    vector<vector<long long>> setup_times(n, vector<long long>(n));
    for (long long i = 0; i < n; ++i) {
        for (long long j = 0; j < n; ++j) {
            long long setup_time = s_dist(generator);
            setup_times[i][j] = (i == j) ? 0 : setup_time;
        }
    }
    long long U_SHARE = U_TARGET / (n / 2);
    long long U_LEFTOVER = U_TARGET - U_SHARE * (n / 2);

    vector<long long> due_dates(n);

    vector<long long> cur_times(NUM_MACHINES, 0);
    vector<long long> end_times(n, 0);
    for (long long i = 0; i < n; ++i) {
        long long prev_av_time = 0;
        for (long long m = 0; m < NUM_MACHINES; ++m) {
            if(i == 0 || m == NUM_MACHINES - 1) {
                processing_times[i][m] = p_dist(generator);
            }else {
                long long base = processing_times[i-1][m+1];
                processing_times[i][m] = max<long long>(base + p_dev_dist(generator), 1);
            }
            if(i != 0)
                cur_times[m] += setup_times[i-1][i];
            cur_times[m] = std::max(prev_av_time, cur_times[m]);
            cur_times[m] += processing_times[i][m];
            prev_av_time = cur_times[m];
        }
        end_times[i] = cur_times.back();
        due_dates[i] = end_times[i];
        if (i % 2 == 1) {
            due_dates[i] = end_times[i] - U_SHARE;
            due_dates[i-1] = end_times[i] + 1;
            setup_times[i][i-1] = max(setup_times[i-1][i], setup_times[i][i-1]) + U_SHARE;
        }
    }
    if ((n/2) % 2 == 1) {
        due_dates[n/2] -= U_LEFTOVER;
    }else {
        due_dates[n/2-1] -= U_LEFTOVER;
    }

    uniform_int_distribution<long long> p_index(0, n/2);
    for(long long i = 0; i < n/50; i++) {
        long long a = p_index(generator);
        long long b = p_index(generator);
        if(a*2 >= n || b*2 >= n)
            continue;
        due_dates[a*2] = std::max(due_dates[a*2], due_dates[b*2]);
    }

    vector<long long> p_indices(n);
    vector<long long> p_indices_rev(n);
    for (long long i = 0; i < n; ++i) {
        p_indices[i] = i;
    }
    shuffle(p_indices.begin(), p_indices.end(), generator);
    for (long long i = 0; i < n; ++i) {
        p_indices_rev[p_indices[i]] = i;
    }
    vector<vector<long long>> processing_times_shuffled(n, vector<long long>(NUM_MACHINES));
    vector<long long> due_dates_shuffled(n);

    for (long long i = 0; i < n; ++i) {
        long long original_index = p_indices[i];
        processing_times_shuffled[i] = processing_times[original_index];
        due_dates_shuffled[i] = due_dates[original_index];
    }
    
    outfile.seekp(0);
    outfile << n << endl;
    
    for (long long i = 0; i < n; ++i) {
        for (long long m = 0; m < NUM_MACHINES; ++m) {
            outfile << processing_times_shuffled[i][m] << " ";
        }
        outfile << due_dates_shuffled[i] << endl; 
    }

    for (long long i = 0; i < n; ++i) {
        long long original_i = p_indices[i];
        
        for (long long j = 0; j < n; ++j) {
            long long original_j = p_indices[j];
            long long mapped_setup_time = setup_times[original_i][original_j]; 
            outfile << mapped_setup_time << (j == n - 1 ? "" : " ");
        }
        outfile << endl;
    }

    outfile.close();
    if (!outfile_sol.is_open()) {
        throw runtime_error("Nie można utworzyć pliku opt");
    }
    outfile_sol.seekp(0);
    outfile_sol << U_TARGET << endl;
    for (long long i = 0; i < n; ++i) {
        outfile_sol << p_indices_rev[i]+1 << " ";
    }
    outfile_sol.close();
}

int main(int argc, char* argv[]) {
    if (argc == 3) {
        string instance_file = argv[1];
        string result_file = argv[2];

        try {
            ProblemData data;
            
            cout << "--- TRYB: WERYFIKACJA ---" << endl;
            read_instance_file(instance_file, data);

            long long expected_tardiness;
            vector<long long> schedule;
            read_result_file(result_file, data.n, expected_tardiness, schedule);

            long long actual_tardiness = evaluate_schedule(data, schedule);

            cout << "\n--- WYNIK WERYFIKACJI ---" << endl;
            cout << "Oczekiwane D: " << expected_tardiness << endl;
            cout << "Obliczone D:   " << actual_tardiness << endl;

            if (actual_tardiness == expected_tardiness) {
                cout << "\nSUCCESS!" << endl;
                return 0;
            } else {
                cout << "\nFAILURE!" << endl;
                cout << "Wartość D w pliku wynikowym (" << expected_tardiness 
                     << ") NIE ZGADZA SIĘ z obliczoną wartością (" << actual_tardiness << ")." << endl;
                return 2;
            }

        } catch (const runtime_error& e) {
            cerr << "\nERR: " << e.what() << endl;
            return 1;
        }
    } else if (argc == 4 && string(argv[1]) == "generate") {
        try {
            cout << "--- TRYB: GENERATOR INSTANCJI ---" << endl;
            string output_file = argv[2];
            size_t n = stoi(argv[3]);
            
            generate_instance(output_file, n);
            return 0;

        } catch (const invalid_argument& e) {
            cerr << "\nARGUMENTY: Wymagana liczba zadań (n) musi być liczbą całkowitą." << endl;
            return 1;
        } catch (const runtime_error& e) {
            cerr << "\nERR: " << e.what() << endl;
            return 1;
        }
    } else {
        cerr << "Użycie: " << endl;
        cerr << "  1. Tryb WERYFIKACJI: " << argv[0] << " <plik_instancji_in> <plik_wynikowy_out>" << endl;
        cerr << "  2. Tryb GENERATORA: " << argv[0] << " generate <plik_wyjściowy.in> <liczba_zadań_n>" << endl;
        return 1;
    }
}
