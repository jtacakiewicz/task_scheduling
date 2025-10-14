#include <cassert>
#include <iomanip>
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
    int sol_id;
    int cmax;
    std::pair<int, int> ids;
    std::pair<int, int> batch_ids;
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

void perform_swap(Solution_t& solution, int id1, int id2, int b1, int b2) {
    auto &vec1 = solution[b1];
    auto &vec2 = solution[b2];

    auto itr1 = std::find(vec1.begin(), vec1.end(), id1);
    auto itr2 = std::find(vec2.begin(), vec2.end(), id2);
    assert(itr1 != vec1.end());
    assert(itr2 != vec2.end());
    auto tmp = *itr1;
    *itr1 = *itr2;
    *itr2 = tmp;
}

bool good_swap(int id1, int id2, const std::vector<std::vector<bool>> &tabu_list) {
    if (id1 == id2)
        return false;
    if (tabu_list[id1].size() > id2 && tabu_list[id1][id2])
        return false;
    else
        return true;
    return true;
}

void add_to_tabu(std::vector<std::vector<bool>> &tabu, int id1, int id2, int n) {
    if (tabu[id1].empty())
        tabu[id1].resize(n, false);
    tabu[id1][id2] = true;
}

Solution_t beam_search(int setuptime, int batch_size, Solution_t batches,
                                const std::vector<JobInfo> &job_info,
                                int max_num_beams=7, int max_num_random=5, double time_limit=5)
{
    double tune_time = 0.01;
    double iters_per_second = -1;
    auto start = std::chrono::steady_clock::now();
    int n = job_info.size();
    int iter_size = n;
    int global_tdur = n;
    int c = calc_cmax(setuptime, batch_size, batches, job_info);

    std::vector<int> global_tabu_list(n*n, -0xffffff);
    std::vector<Solution_t> old_winners;
    std::vector<Solution_t> winners = {batches};
    std::vector<Swap> last_swaps;
    int global_iter = 0;
    int best_cmax_so_far = 0xffffff;
    Solution_t best_winner;

    while (true) {
        auto now = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double>(now - start).count();
        if (elapsed > tune_time && iters_per_second == -1) {
            iters_per_second = global_iter / elapsed;
            std::cerr << "Search iterations per second: " << iters_per_second << std::endl;
        }
        if (elapsed > time_limit)
            break;
        global_iter++;

        float t = 0;
        if(iters_per_second != -1)
            t = global_iter / float(iters_per_second*time_limit);
        t = sqrt(t);
        int num_beams = t * max_num_beams;
        int num_rand = (1.f - t) * max_num_random;

        std::vector<Swap> swaps;


        for (int sol_id=0; sol_id < winners.size(); sol_id++) {
            auto& solution = winners[sol_id];
            std::vector<int> batch_id(n, -1);
            for (int i = 0; i < (int)solution.size(); i++)
                for (int job : solution[i]){ 
                    batch_id[job] = i;
                }

            std::vector<std::vector<bool>> local_tabu_list(n);

            for (int local_iter = 0; local_iter < iter_size; local_iter++) {
                int id1 = 0, id2 = 0;
                auto b1 = batch_id[id1];
                auto b2 = batch_id[id2];
                while (batch_id[id1] == batch_id[id2]
                       || !good_swap(id1, id2, local_tabu_list)
                       || global_tabu_list[id1*batch_size+b2]>global_iter - global_tdur 
                       || global_tabu_list[id2*batch_size+b1]>global_iter - global_tdur)
                {
                    int i1 = rand() % n;
                    int i2 = rand() % n;
                    id1 = std::min(i1, i2);
                    id2 = std::max(i1, i2);
                    b1 = batch_id[id1];
                    b2 = batch_id[id2];
                }
                add_to_tabu(local_tabu_list, id1, id2, n);

                perform_swap(solution, id1, id2, b1, b2);
                int cmax = calc_cmax(setuptime, batch_size, solution, job_info);
                perform_swap(solution, id1, id2, b2, b1);
                swaps.push_back({sol_id, cmax, {id1, id2}, {b1, b2}});
            }
        }

        // swaps.insert(swaps.end(), last_swaps.begin(), last_swaps.end());
        sort(swaps.begin(), swaps.end(), [](auto &a, auto &b){
            return a.cmax < b.cmax;
        });

        last_swaps.clear();
        old_winners.clear();
        std::swap(old_winners, winners);

        int max_win_id = std::min((int)swaps.size(), num_beams);
        std::vector<int> win_idx(max_win_id);
        iota(win_idx.begin(), win_idx.end(), 0);

        std::vector<int> rand_idx;
        for (int i = 0; i < num_rand && max_win_id < (int)swaps.size(); i++)
            rand_idx.push_back(rand() % ((int)swaps.size() - max_win_id) + max_win_id);

        std::vector<int> chosen = win_idx;
        chosen.insert(chosen.end(), rand_idx.begin(), rand_idx.end());

        for (int i : chosen) {
            const auto& swap = swaps[i];
            last_swaps.push_back(swap);
            Solution_t sol = old_winners[swap.sol_id];
            perform_swap(sol, swap.ids.first, swap.ids.second, swap.batch_ids.first, swap.batch_ids.second);

            int id1 = swap.ids.first;
            int id2 = swap.ids.second;
            int b1 = swap.batch_ids.first;
            int b2 = swap.batch_ids.second;
            global_tabu_list[id1*batch_size+b1] = global_iter;
            global_tabu_list[id2*batch_size+b2] = global_iter;

            if (swap.cmax < best_cmax_so_far) {
                std::cerr << "New best: " << swap.cmax << std::endl;
                best_cmax_so_far = swap.cmax;
                best_winner = sol;
            }
            winners.push_back(sol);
        }
    }

    return std::move(best_winner);
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

    float avg_task = 0;
    for(auto j : job_info) {
        avg_task += j.duration;
    }
    avg_task /= job_info.size();

    auto solution = solve_greedy(setup, batch_s, job_info);
    int cgreedy = calc_cmax(setup, batch_s, solution, job_info);
    std::cerr << "greedy solution: " << cgreedy  << std::endl;
    solution = beam_search(setup, batch_s, solution, job_info);
    int ctabu = calc_cmax(setup, batch_s, solution, job_info);
    std::cerr << "search solution: " << ctabu << std::endl;
    std::cerr << "search improvement: " << std::setprecision(2) << (cgreedy-ctabu)/avg_task*100 <<"% of a task\n";

    std::cout << "0" << std::endl; // cmax
    std::cout << solution.size() << std::endl; // number of batches
    for (auto &b : solution) {
        for (int j : b)
            std::cout << j + 1 << " ";
        std::cout << std::endl;
    }

    return 0;
}

