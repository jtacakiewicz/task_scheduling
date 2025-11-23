#include <iomanip>
#include <iterator>
#include <random>
#include <unordered_set>
#include <vector>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <math.h>
#include <iostream>
#include <chrono>
#include <climits>
#include "thread_pool.hpp"

typedef unsigned long long metric_t;
struct JobInfo {
    metric_t ready;
    metric_t duration;
    metric_t due;
};
void handle_alarm() {
    std::cout << "alarmed";
}
std::vector<JobInfo> load_in_file(std::istream &inp) {
    metric_t n;
    inp >> n;
    if (n <= 0)
        throw std::invalid_argument("Invalid input parameters");

    std::vector<JobInfo> job_info(n);
    for (int i = 0; i < n; i++) {
        metric_t duration, ready_time, due_time;
        inp >> duration >> ready_time >> due_time;
        if (duration < 0 || ready_time < 0 || due_time < 0)
            throw std::invalid_argument("Job properties cannot be less than zero");
        job_info[i] = {ready_time, duration, due_time};
    }
    return job_info;
}


int main(int argc, char **argv) {
    auto start = std::chrono::steady_clock::now();
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <input_file> <output_file> <time_limit_sec>" << std::endl;
        return 1;
    }

    std::ifstream fin(argv[1]);
    std::ofstream fout(argv[2]);
    double sec_limit = std::stod(argv[3]);
    if (!fin) {
        std::cerr << "Error opening file: " << argv[1] << std::endl;
        return 1;
    }

    int setup, batch_s;
    std::vector<JobInfo> job_info;
    tie(setup, batch_s, job_info) = load_in_file(fin);

    float avg_task = 0;
    for(auto j : job_info) {
        avg_task += j.duration;
    }
    avg_task /= job_info.size();

    auto solution = solve_greedy(setup, batch_s, job_info);
    metric_t cgreedy = calc_cmax(setup, batch_s, solution, job_info);

    auto now = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(now - start).count();
    std::cerr << "Time spend up to beam search: " << elapsed << std::endl;
    sec_limit -= elapsed;
    sec_limit *= 0.95;

    std::cerr << "greedy solution: " << cgreedy  << std::endl;
    solution = beam_search(setup, batch_s, solution, job_info, sec_limit);
    metric_t ctabu = calc_cmax(setup, batch_s, solution, job_info);
    std::cerr << "search solution: " << ctabu << std::endl;
    std::cerr << "search improvement: " << std::setprecision(2) << (cgreedy-ctabu)/avg_task*100 <<"% of a task\n";

    fout << ctabu << std::endl; // cmax
    fout << solution.size() << std::endl; // number of batches
    for (auto &b : solution) {
        for (int j : b)
            fout << j + 1 << " ";
        fout << std::endl;
    }

    return 0;
}

