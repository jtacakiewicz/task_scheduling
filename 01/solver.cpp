#include <vector>
#include <optional>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <math.h>
#include <iostream>
#include <chrono>
#include <climits>

typedef std::vector<std::vector<int>> Solution_t;

struct JobInfo {
    int ready;
    int duration;
};
struct Swap {
    int cmax;
    Solution_t solution;
    int job1;
    int job2;
};

int finish_time(const JobInfo &job) {
    return job.ready + job.duration;
}

std::tuple<int,int,std::vector<JobInfo>> check_in_file(std::istream &inp) {
    int n, setuptime, batch_size;
    inp >> n >> setuptime >> batch_size;
    fprintf(stderr, "n: %d; setup: %d; batch: %d\n", n, setuptime, batch_size);

    if (n <= 0 || setuptime < 0 || batch_size <= 0)
        throw std::invalid_argument("Invalid input parameters");

    std::vector<JobInfo> job_info(n);
    for (int i = 0; i < n; i++) {
        int processing_time, ready_time;
        inp >> processing_time >> ready_time;
        if (processing_time < 0)
            throw std::invalid_argument("Processing time cannot be less than zero");
        if (ready_time < 0)
            throw std::invalid_argument("Ready time cannot be less than zero");
        job_info[i] = {ready_time, processing_time};
    }
    return {setuptime, batch_size, job_info};
}

Solution_t solve_greedy(int setuptime, int batch_size, const std::vector<JobInfo> &job_info) {
    std::vector<int> jobs(job_info.size());
    iota(jobs.begin(), jobs.end(), 0);
    sort(jobs.begin(), jobs.end(), [&](int a, int b){
        return job_info[a].ready > job_info[b].ready;
    });

    int time = job_info[jobs[jobs.size() - batch_size]].ready;
    std::cerr << "Time as the start: " << time << std::endl;
    Solution_t batches;

    while (!jobs.empty()) {
        if (batches.empty() || !batches.back().empty())
            batches.push_back({});
        int left = (int)jobs.size() - 1;
        while (left >= 0 && job_info[jobs[left]].ready <= time)
            left--;
        std::vector<int> potential(jobs.begin() + left + 1, jobs.end());
        sort(potential.begin(), potential.end(), [&](int a, int b){
            return job_info[a].duration < job_info[b].duration;
        });

        int new_time = time;
        int cur_batch_size = 0;
        for (int job : potential) {
            if(cur_batch_size == batch_size)
                break;
            if (job_info[job].ready > time)
                break;
            new_time = std::max(time + job_info[job].duration, new_time);
            batches.back().push_back(job);
            jobs.erase(remove(jobs.begin(), jobs.end(), job), jobs.end());
            cur_batch_size += 1;
        }
        time = new_time + setuptime;
    }

    if (!batches.empty() && batches.back().empty())
        batches.pop_back();

    return batches;
}

int calc_cmax(int setuptime, int batch_size, const Solution_t &batches, const std::vector<JobInfo> &job_info) {
    int cmax = 0;
    for (auto &batch : batches) {
        int ready = 0, duration = 0;
        for (int job : batch) {
            ready = std::max(ready, job_info[job].ready);
            duration = std::max(duration, job_info[job].duration);
        }
        cmax = std::max(ready, cmax);
        cmax += duration;
        cmax += setuptime;
    }
    cmax -= setuptime;
    return cmax;
}

Solution_t perform_swap(Solution_t solution, const std::vector<int> &batch_id, int id1, int id2) {
    int b1 = batch_id[id1];
    int b2 = batch_id[id2];

    auto &vec1 = solution[b1];
    auto &vec2 = solution[b2];

    vec1.erase(remove(vec1.begin(), vec1.end(), id1), vec1.end());
    vec2.erase(remove(vec2.begin(), vec2.end(), id2), vec2.end());
    vec1.push_back(id2);
    vec2.push_back(id1);
    return solution;
}

bool good_swap(int i1, int i2, const std::vector<std::vector<std::optional<int>>> &tabu_list, int min_iter=-1) {
    if (i1 == i2)
        return false;
    int id1 = std::min(i1, i2);
    int id2 = std::max(i1, i2);
    if (tabu_list[id1].size() > id2 && tabu_list[id1][id2].has_value()) {
        if (tabu_list[id1][id2].value() > min_iter)
            return false;
        else
            return true;
    }
    return true;
}

void add_to_tabu(std::vector<std::vector<std::optional<int>>> &tabu, int id1, int id2, int n, int iter=1) {
    if (tabu[id1].empty())
        tabu[id1].resize(n);
    tabu[id1][id2] = iter;
}

Solution_t beam_search(int setuptime, int batch_size, Solution_t batches,
                                const std::vector<JobInfo> &job_info,
                                int num_beams=7, int num_rand=2, double time_limit=5)
{
    double tune_time = 0.1;
    double iters_per_second = -1;
    auto start = std::chrono::steady_clock::now();
    int n = job_info.size();
    int iter_size = n;
    int global_tdur = 100;
    int c = calc_cmax(setuptime, batch_size, batches, job_info);

    std::vector<std::vector<std::optional<int>>> global_tabu_list(n);
    std::vector<Solution_t> winners = {batches};
    std::vector<Swap> last_swaps;
    int global_iter = 0;
    int best_cmax_so_far = INT_MAX;
    Solution_t best_winner;

    while (true) {
        auto now = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double>(now - start).count();
        if (elapsed > tune_time && iters_per_second == -1) {
            iters_per_second = global_iter / elapsed;
            std::cerr << "Search iterations per second: " << iters_per_second << std::endl;
            global_tdur = global_iter * 0.1;
        }
        if (elapsed > time_limit)
            break;
        global_iter++;

        std::vector<Swap> swaps;

        for (auto &solution : winners) {
            std::vector<int> batch_id(n, -1);
            for (int i = 0; i < (int)solution.size(); i++)
                for (int job : solution[i])
                    batch_id[job] = i;

            std::vector<std::vector<std::optional<int>>> local_tabu_list(n);

            for (int local_iter = 0; local_iter < iter_size; local_iter++) {
                int i1 = rand() % n;
                int i2 = i1;
                while (batch_id[i1] == batch_id[i2]
                       || !good_swap(i1, i2, local_tabu_list)
                       || !good_swap(i1, i2, global_tabu_list, global_iter - global_tdur))
                {
                    i2 = rand() % n;
                }

                int id1 = std::min(i1, i2);
                int id2 = std::max(i1, i2);
                add_to_tabu(local_tabu_list, id1, id2, n);

                auto temp = perform_swap(solution, batch_id, id1, id2);
                int cmax = calc_cmax(setuptime, batch_size, temp, job_info);
                swaps.push_back({cmax, temp, id1, id2});
            }
        }

        swaps.insert(swaps.end(), last_swaps.begin(), last_swaps.end());
        sort(swaps.begin(), swaps.end(), [](auto &a, auto &b){
            return a.cmax < b.cmax;
        });

        last_swaps.clear();
        winners.clear();

        int max_win_id = std::min((int)swaps.size(), num_beams);
        std::vector<int> win_idx(max_win_id);
        iota(win_idx.begin(), win_idx.end(), 0);

        std::vector<int> rand_idx;
        for (int i = 0; i < num_rand && max_win_id < (int)swaps.size(); i++)
            rand_idx.push_back(rand() % ((int)swaps.size() - max_win_id) + max_win_id);

        std::vector<int> chosen = win_idx;
        chosen.insert(chosen.end(), rand_idx.begin(), rand_idx.end());

        for (int i : chosen) {
            last_swaps.push_back(swaps[i]);
            winners.push_back(swaps[i].solution);
            int id1 = swaps[i].job1;
            int id2 = swaps[i].job2;
            add_to_tabu(global_tabu_list, id1, id2, n, global_iter);
        }

        if (swaps[0].cmax < best_cmax_so_far) {
            best_cmax_so_far = swaps[0].cmax;
            best_winner = swaps[0].solution;
        }

        if(iters_per_second != -1 && global_iter % int(iters_per_second) == 0)
            std::cerr << "solution of iter " << global_iter << ": " << swaps[0].cmax << std::endl;
    }

    return best_winner;
}

int main(int argc, char **argv) {
    srand(42);
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
        return 1;
    }

    std::ifstream fin(argv[1]);
    if (!fin) {
        std::cerr << "Error opening file." << std::endl;
        return 1;
    }

    int setup, batch_s;
    std::vector<JobInfo> job_info;
    tie(setup, batch_s, job_info) = check_in_file(fin);

    auto solution = solve_greedy(setup, batch_s, job_info);
    std::cerr << "greedy solution: " << calc_cmax(setup, batch_s, solution, job_info) << std::endl;
    solution = beam_search(setup, batch_s, solution, job_info);

    std::cout << "0" << std::endl; // cmax
    std::cout << solution.size() << std::endl; // number of batches
    for (auto &b : solution) {
        for (int j : b)
            std::cout << j + 1 << " ";
        std::cout << std::endl;
    }

    std::cerr << calc_cmax(setup, batch_s, solution, job_info) << std::endl;
    return 0;
}

