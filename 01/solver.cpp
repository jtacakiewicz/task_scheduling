#include <cassert>
#include <cstdlib>
#include <cstring>
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

typedef std::vector<std::vector<int>> Solution_t;

typedef unsigned long long metric_t;
struct JobInfo {
    metric_t ready;
    metric_t duration;
};
struct Swap {
    int sol_id;
    metric_t cmax;
    std::pair<int, int> ids;
    std::pair<int, int> batch_ids;
};

int finish_time(const JobInfo &job) {
    return job.ready + job.duration;
}

#define MAX_WINNERS 30
#define MAXN 500
#define INVALID_ID (MAXN*2)

std::tuple<int,int,std::vector<JobInfo>> load_in_file(std::istream &inp) {
    metric_t n, setuptime, batch_size;
    inp >> n >> setuptime >> batch_size;
    fprintf(stderr, "n: %llu; setup: %llu; batch: %llu\n", n, setuptime, batch_size);

    if (n <= 0 || setuptime < 0 || batch_size <= 0)
        throw std::invalid_argument("Invalid input parameters");

    std::vector<JobInfo> job_info(n);
    for (int i = 0; i < n; i++) {
        metric_t processing_time, ready_time;
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

    auto time = job_info[jobs.back()].ready;
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

        auto new_time = time;
        auto cur_batch_size = 0;
        for (int job : potential) {
            if(cur_batch_size == batch_size)
                break;
            if (job_info[job].ready > time) {
                break;
            }
            new_time = std::max(time + job_info[job].duration, new_time);
            batches.back().push_back(job);
            jobs.erase(remove(jobs.begin(), jobs.end(), job), jobs.end());
            cur_batch_size += 1;
        }
        time = std::max(new_time + setuptime, job_info[*(jobs.begin() + left)].ready);
    }

    if (!batches.empty() && batches.back().empty())
        batches.pop_back();

    return batches;
}

metric_t calc_cmax(int setuptime, int batch_size, const Solution_t &batches, const std::vector<JobInfo> &job_info) {
    metric_t cmax = 0;
    for (auto &batch : batches) {
        metric_t ready = 0, duration = 0;
        for (int job : batch) {
            ready = std::max(ready, job_info[job].ready);
            duration = std::max(duration, job_info[job].duration);
        }
        cmax = std::max(ready, cmax);
        cmax += duration;
        if(batch.size() != 0) {
            cmax += setuptime;
        }
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
void perform_move(Solution_t& solution, int id, int b1, int b2) {
    auto& batchA = solution[b1];
    auto& batchB = solution[b2];
    auto itr = std::find(batchA.begin(), batchA.end(), id);
    assert(itr != batchA.end());
    batchA.erase(itr);
    batchB.push_back(id);
}

//the higher the lambda the more close to the top the results will be
std::vector<int> pick_biased_indices(size_t n, int pick_winners, int pick_rand, std::mt19937& gen, double lmbda=0.05) {
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
typedef int16_t tabu_list[MAXN * MAXN];
struct CensoredSolution {
    Solution_t conf;
    tabu_list tabu_ind;
    tabu_list tabu_batch;
};
CensoredSolution old_winners[MAX_WINNERS];
CensoredSolution winners[MAX_WINNERS];
ThreadPool thread_pool;
std::vector<Swap> search_moves(CensoredSolution& solution, int sol_id, int n, int global_iter, int tabu_ind, int tabu_batch, int setuptime, int batch_size, const std::vector<JobInfo> &job_info, int check_per_conf, std::mt19937::result_type seed) {
    std::mt19937 gen(seed);
    std::vector<Swap> local_swaps;
    std::vector<int> batch_id(n, -1);
    for (int i = 0; i < (int)solution.conf.size(); i++)
        for (int job : solution.conf[i])
            batch_id[job] = i;
    std::vector<int> free_batches;
    for (int i = 0; i < (int)solution.conf.size(); i++)
        if(solution.conf[i].size() < batch_size)
            free_batches.push_back(i);
    check_per_conf *= free_batches.size();
        
    int tries = 0;
    const int max_tries = check_per_conf * 100;
    std::uniform_int_distribution<> id_rand(0, n - 1);
    std::uniform_int_distribution<> batch_rand(0, free_batches.size() - 1);
    std::unordered_set<int> performed_moves;
    while(local_swaps.size() < check_per_conf && tries < max_tries && !free_batches.empty()) {
        ++tries;
        int id = id_rand(gen);
        auto b = batch_id[id];
        int b_other = batch_rand(gen);
        b_other = free_batches[b_other];
        auto hash = (id * 7631563) ^ (b_other*5928253);
        if(b == b_other) {
            continue;
        }
        if(performed_moves.find(hash) != performed_moves.end()) {
            continue;
        }
        performed_moves.insert(hash);
        //perform swap
        perform_move(solution.conf, id, b, b_other);
        auto cmax = calc_cmax(setuptime, batch_size, solution.conf, job_info);
        //roll back
        perform_move(solution.conf, id, b_other, b);
        local_swaps.push_back({sol_id, cmax, {id, INVALID_ID}, {b, b_other}});
    }
    return local_swaps;
}
std::vector<Swap> search_swaps(CensoredSolution& solution, int sol_id, int n, int global_iter, int tabu_ind, int tabu_batch, int setuptime, int batch_size, const std::vector<JobInfo> &job_info, int check_per_conf, std::mt19937::result_type seed) {
    std::mt19937 gen(seed);
    std::vector<Swap> local_swaps;
    std::vector<int> batch_id(n, -1);
    for (int i = 0; i < (int)solution.conf.size(); i++)
        for (int job : solution.conf[i])
            batch_id[job] = i;

    int tries = 0;
    const int max_tries = check_per_conf * 100;
    int found = 0;
    std::uniform_int_distribution<> dist(0, n - 1);


    while (found < check_per_conf && tries < max_tries) {
        ++tries;
        int id1 = dist(gen);
        int id2 = dist(gen);
        if (id1 == id2) continue; 
        if (id1 > id2) std::swap(id1, id2);

        auto b1 = batch_id[id1];
        auto b2 = batch_id[id2];
        auto bX = std::max(b1, b2);
        auto bM = std::min(b1, b2);
        auto hash = (id1 * 7631563) ^ (id2*5928253);

        //dont swap within the same batch
        if (b1 == b2)
            continue;
        //respect tabu lists
        if (solution.tabu_ind[id1*MAXN + id2] > global_iter - tabu_ind)
            continue;
        if (solution.tabu_batch[bM*MAXN + bX] > global_iter - tabu_batch)
            continue;

        perform_swap(solution.conf, id1, id2, b1, b2);
        metric_t cmax = calc_cmax(setuptime, batch_size, solution.conf, job_info);
        perform_swap(solution.conf, id1, id2, b2, b1);
        local_swaps.push_back({sol_id, cmax, {id1, id2}, {b1, b2}});
        ++found;
    }
    return local_swaps;
}

Solution_t beam_search(const int setuptime, const int batch_size, Solution_t best_sol,
                                const std::vector<JobInfo> &job_info, double time_limit=5,
                                int pick_winners=3, int pick_rand=5)
{
    std::mt19937 gen(42);
    double tune_time = 1.;
    double iters_per_second = -1;
    auto start = std::chrono::steady_clock::now();
    int n = job_info.size();
    int tabu_dur = n;
    int tabu_batch_dur = best_sol.size();
    int check_per_conf = n * 2;
    int move_checks = check_per_conf * 0.1;

    int old_winners_size = 0;
    int winners_size = 1;
    winners[0].conf = best_sol;
    for(int i = 0; i < MAX_WINNERS; i++) {
        std::fill(std::begin(old_winners[i].tabu_ind), std::end(old_winners[i].tabu_ind), INT16_MIN);
        std::fill(std::begin(winners[i].tabu_ind), std::end(winners[i].tabu_ind), INT16_MIN);

        std::fill(std::begin(old_winners[i].tabu_batch), std::end(old_winners[i].tabu_batch), INT16_MIN);
        std::fill(std::begin(winners[i].tabu_batch), std::end(winners[i].tabu_batch), INT16_MIN);
    }

    std::vector<Swap> last_swaps;
    int global_iter = 1;
    metric_t best_cmax_so_far = calc_cmax(setuptime, batch_size, best_sol, job_info);

    while (true) {
        assert(winners_size != 0);
        auto now = std::chrono::steady_clock::now();
        double elapsed = std::chrono::duration<double>(now - start).count();
        double sec_per_iter = elapsed/global_iter;
        if (elapsed > tune_time && iters_per_second == -1) {
            iters_per_second = global_iter / elapsed;
            std::cerr << "Search iterations per second: " << iters_per_second*(pick_winners*pick_rand) << std::endl;
        }
        if (elapsed+sec_per_iter*1.2>time_limit) {
            break;
        }
        global_iter++;

        std::vector<Swap> swaps;

        std::mutex swaps_results_mut;
        for (int sol_id=0; sol_id < winners_size; sol_id++) {
            auto seed = gen();
            thread_pool.addTask([&, sol_id, seed]() {
                auto lswaps = search_swaps(winners[sol_id], sol_id, n, global_iter, tabu_dur, tabu_batch_dur, setuptime, batch_size, job_info, check_per_conf, seed);
                auto lmoves = search_moves(winners[sol_id], sol_id, n, global_iter, tabu_dur, tabu_batch_dur, setuptime, batch_size, job_info, move_checks, seed);
                swaps_results_mut.lock();
                swaps.insert(swaps.end(), lswaps.begin(), lswaps.end());
                swaps.insert(swaps.end(), lmoves.begin(), lmoves.end());
                swaps_results_mut.unlock();
            });
        }
        thread_pool.waitForCompletion();

        if(swaps.size() == 0) {
            std::cerr << "Unexpected behaviour!\n";
            return best_sol;
        }
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
            winners[winners_size] = old_winners[swap.sol_id];
            auto& sol = winners[winners_size];
            if(swap.ids.second != INVALID_ID) {
                perform_swap(sol.conf, swap.ids.first, swap.ids.second, swap.batch_ids.first, swap.batch_ids.second);
                int id1 = swap.ids.first;
                int id2 = swap.ids.second;
                int b1 = swap.batch_ids.first;
                int b2 = swap.batch_ids.second;
                auto bX = std::max(b1, b2);
                auto bM = std::min(b1, b2);
                sol.tabu_ind[id1*MAXN+id2] = global_iter;
                sol.tabu_batch[bM*MAXN+bX] = global_iter;
            }else if (swap.batch_ids.second != INVALID_ID) {
                perform_move(sol.conf, swap.ids.first, swap.batch_ids.first, swap.batch_ids.second);
                if(sol.conf[swap.batch_ids.first].size() == 0) {
                    auto itr = sol.conf.begin() + swap.batch_ids.first;
                    assert(itr->size() == 0);
                    sol.conf.erase(itr);
                }
            }

            if (swap.cmax < best_cmax_so_far) {
                std::cerr << "New best: " << swap.cmax << std::endl;
                best_cmax_so_far = swap.cmax;
                best_sol = sol.conf;
            }
            winners_size += 1;
        }
    }

    return std::move(best_sol);
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

