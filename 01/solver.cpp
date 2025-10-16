#include <cassert>
#include <cstring>
#include <iomanip>
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

#define MAX_WINNERS 30
#define MAXN 500

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

//the higher the lambda the more close to the top the results will be
std::vector<int> pick_biased_indices(size_t n, int pick_winners, int pick_rand, std::mt19937& gen, double lmbda=0.01) {
    int max_win_id = std::min((int)n, pick_winners);

    std::vector<int> chosen;
    chosen.reserve(pick_winners + pick_rand);

    for (int i = 0; i < max_win_id; ++i)
        chosen.push_back(i);
    std::unordered_set<int> used(chosen.begin(), chosen.end());
    std::exponential_distribution<> exp_dist(lmbda);

    while ((int)chosen.size() < pick_winners + pick_rand) {
        double sample = exp_dist(gen);
        int idx = static_cast<int>(sample);

        if (idx < 0 || idx >= (int)n) continue;
        if (used.insert(idx).second)
            chosen.push_back(idx);
    }

    return chosen;
}

Solution_t beam_search(int setuptime, int batch_size, Solution_t best_sol,
                                const std::vector<JobInfo> &job_info,
                                int pick_winners=3, int pick_rand=5, double time_limit=5)
{
    std::mt19937 gen(42);
    double tune_time = 1.;
    double iters_per_second = -1;
    auto start = std::chrono::steady_clock::now();
    int n = job_info.size();
    int global_tdur = best_sol.size() * n;
    int global_dur_batch = best_sol.size();

    int global_tabu_list[MAXN * MAXN]{-MAXN};
    int global_tabu_batch[MAXN * MAXN]{-MAXN};
    for(int i = 0; i < MAXN * MAXN; i++) {
        global_tabu_list[i] = -0xffffff;
        global_tabu_batch[i] = -0xffffff;
    }

    Solution_t old_winners[MAX_WINNERS];
    int old_winners_size = 0;
    Solution_t winners[MAX_WINNERS];
    int winners_size = 1;
    winners[0] = best_sol;

    std::vector<Swap> last_swaps;
    int global_iter = 0;
    int best_cmax_so_far = calc_cmax(setuptime, batch_size, best_sol, job_info);

    while (true) {
        assert(winners_size != 0);
        auto now = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double>(now - start).count();
        if (elapsed > tune_time && iters_per_second == -1) {
            iters_per_second = global_iter / elapsed;
            std::cerr << "Search iterations per second: " << iters_per_second << std::endl;
        }
        if (elapsed > time_limit)
            break;
        global_iter++;

        std::vector<Swap> swaps;

        for (int sol_id=0; sol_id < winners_size; sol_id++) {
            auto& solution = winners[sol_id];
            std::vector<int> batch_id(n, -1);
            for (int i = 0; i < (int)solution.size(); i++)
                for (int job : solution[i])
                    batch_id[job] = i;

            for(int id1 = 0; id1 < n; id1++) {
                for(int id2 = id1+1; id2 < n; id2++) {
                    auto b1 = batch_id[id1];
                    auto b2 = batch_id[id2];
                    auto bX = std::max(b1, b2);
                    auto bM = std::min(b1, b2);
                    if(batch_id[id1] == batch_id[id2])
                        continue;
                    if(global_tabu_list[id1*n+id2]>global_iter - global_tdur )
                        continue;
                    if(global_tabu_batch[bM*n+bX]>global_iter - global_dur_batch )
                        continue;
                    perform_swap(solution, id1, id2, b1, b2);
                    int cmax = calc_cmax(setuptime, batch_size, solution, job_info);
                    perform_swap(solution, id1, id2, b2, b1);
                    swaps.push_back({sol_id, cmax, {id1, id2}, {b1, b2}});
                }
            }
        }

        // for(auto& l : last_swaps) {
        //     winners[winners_size] = old_winners[l.sol_id];
        //     l.sol_id = winners_size;
        //     winners_size++;
        // }

        // swaps.insert(swaps.end(), last_swaps.begin(), last_swaps.end());
        sort(swaps.begin(), swaps.end(), [](auto &a, auto &b){
            return a.cmax < b.cmax;
        });

        last_swaps.clear();
        old_winners_size = 0;
        std::swap(old_winners, winners);
        std::swap(old_winners_size, winners_size);

        auto chosen = pick_biased_indices(swaps.size(), pick_winners, pick_rand, gen);

        for (int i : chosen) {
            const auto& swap = swaps[i];
            last_swaps.push_back(swap);
            assert(old_winners_size > swap.sol_id);
            Solution_t sol = old_winners[swap.sol_id];
            perform_swap(sol, swap.ids.first, swap.ids.second, swap.batch_ids.first, swap.batch_ids.second);

            int id1 = swap.ids.first;
            int id2 = swap.ids.second;
            int b1 = swap.batch_ids.first;
            int b2 = swap.batch_ids.second;
            auto bX = std::max(b1, b2);
            auto bM = std::min(b1, b2);
            global_tabu_list[id1*n+id2] = global_iter;
            global_tabu_batch[bM*n+bX] = global_iter;

            if (swap.cmax < best_cmax_so_far) {
                std::cerr << "New best: " << swap.cmax << std::endl;
                best_cmax_so_far = swap.cmax;
                best_sol = sol;
            }
            winners[winners_size] = sol;
            winners_size += 1;
        }
    }

    return std::move(best_sol);
}

int main(int argc, char **argv) {
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

