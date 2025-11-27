#include <cassert>
#include <array>
#include <iterator>
#include <numeric>
#include <random>
#include <set>
#include <vector>
#include <fstream>
#include <math.h>
#include <iostream>
#include <chrono>
#include <climits>
#include <stdio.h>
#include <unistd.h>
#include <signal.h>
#include <stdlib.h>
#include <algorithm>
#include "tree.hpp"

typedef unsigned long long metric_t;
typedef unsigned short JobID_t;
#define NUM_MACHINES 5
#define MAX_JOBS 500
typedef std::array<std::vector<JobID_t>, NUM_MACHINES> Config_t;

std::mt19937 rng(std::random_device{}());
struct JobInfo {
    metric_t ready;
    metric_t duration;
    metric_t due;
};
void alarm_handler(int signum) {
    (void)signum;
    exit(0);
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

std::tuple<int, Config_t, std::vector<JobID_t>> solve_greedy(const std::vector<JobInfo>& job_info, float coef = 1.f) {
    std::vector<JobID_t> jobs(job_info.size());
    std::iota(jobs.begin(), jobs.end(), 0);
    std::sort(jobs.begin(), jobs.end(), [&](JobID_t a, JobID_t b) {
        return job_info[a].ready < job_info[b].ready;
    });

    double avg_duration = 0;
    for(auto i : job_info) avg_duration += i.duration;
    avg_duration /= job_info.size();


    std::array<metric_t, NUM_MACHINES> machine_times;
    Config_t result;
    std::fill(machine_times.begin(), machine_times.end(), 0);

    std::vector<JobID_t> to_readd;
    int U = 0;
    for(auto job_id : jobs) {
        if(job_info[job_id].duration > avg_duration * coef) {
            to_readd.push_back(job_id);
            continue;
        }
        auto& job = job_info[job_id];
        size_t best_machine = -1;
        metric_t best_start = -1;
        for (int i = 0; i < NUM_MACHINES; i++) {
            metric_t start = std::max(job.ready, machine_times[i]);
            if(best_start == -1 || start < best_start) {
                best_start = start;
                best_machine = i;
            }
        }
        assert(best_start != -1);
        metric_t start = best_start;
        metric_t finish = start + job.duration;
        machine_times[best_machine] = finish;
        result[best_machine].push_back(job_id);
        if(finish > job.due)
            U += 1;
    }
    return std::make_tuple(U, result, to_readd);
}
struct Range {
    metric_t l, r;
    long long value;
};

// get all ranges where for the entire range the amount of ready jobs is not less than threshold
std::vector<Range> get_ranges(const std::vector<JobInfo>& job_info, const std::vector<JobID_t>& jobs, size_t threshold_min, float operating_coef=1.0, size_t threshold_max = 0xffffff) {
    std::vector<std::pair<metric_t, int>> events; 
    events.reserve(jobs.size() * 2);

    for (auto &j : jobs) {
        auto& job = job_info[j];
        events.push_back({job.ready, +1});
        events.push_back({job.due - job.duration*operating_coef, -1});  // last possible moment
    }
    sort(events.begin(), events.end());

    std::vector<Range> result;
    long long active = 0;

    auto isValid = [&](long long ctr) {
        return !(ctr < threshold_min || ctr > threshold_max);
    };

    metric_t seg_start = 0;
    for (size_t i = 0; i < events.size(); ++i) {
        metric_t x = events[i].first;
        int delta = events[i].second;

        long long prev_active = active;
        active += delta;
        while(events[i].first == x && i+1 < events.size()) {
            i++;
            active += events[i].second;
        }

        if (isValid(active) && !isValid(prev_active)) {
            seg_start = x;
        }

        if (isValid(prev_active) && !isValid(active)) {
            metric_t seg_end = x;
            result.push_back({seg_start, seg_end, prev_active});
        }
    }
    assert(active == 0);

    return result;
}
// get most number of ranges based on simulated annealing
std::vector<Range> search_ranges(const std::vector<JobInfo>& job_info, const std::vector<JobID_t>& jobs, int max_iters=10000) {
    metric_t max_due = 0;
    for(auto j : jobs) {
        max_due = std::max(max_due, job_info[j].due - job_info[j].duration);
    }
    float t_init = 0.5f;
    float t_min = 0.0f;
    float t_max = 1.0f;
    std::uniform_real_distribution<float> dist(-0.25f, 0.25f); // step size
    std::uniform_real_distribution<float> prob(0.0f, 1.0f);

    float current_t = t_init;
    long long current_score = get_ranges(job_info, jobs, (max_due * current_t)).size();

    float temperature = 1.0f;
    float cooling_rate = 0.999f;

    float best_t = current_t;
    long long best_score = current_score;

    int iter = 0;
    for (; iter < max_iters; ++iter) {
        // propose a new t by small random step
        float new_t = current_t + dist(rng) * (0.1f + temperature);
        new_t = std::clamp(new_t, t_min, t_max);

        long long new_score = get_ranges(job_info, jobs, (max_due * new_t)).size();
        long long delta = new_score - current_score;

        // accept if better, or with probability exp(delta/T)
        if (delta >= 0 || prob(rng) < std::exp(delta / temperature)) {
            current_t = new_t;
            current_score = new_score;

            if (new_score > best_score) {
                best_t = new_t;
                best_score = new_score;
            }
        }

        temperature *= cooling_rate; // cool down
        if (temperature < 1e-6) break;
    }
    return get_ranges(job_info, jobs, (max_due * best_t));
}
std::vector<float> getPrio(const std::vector<JobInfo>& job_info) {
    std::vector<JobID_t> jobs(job_info.size());
    std::iota(jobs.begin(), jobs.end(), 0);
    std::vector<float> priority(job_info.size(), 0.f);
    auto ranges = search_ranges(job_info, jobs);
    for(auto r : ranges) {
        for(auto j : jobs) {
            auto& job = job_info[j];
            if(!(job.due - job.duration < r.l || r.r < job.ready)) {
                priority[j] += 1;
            }
        }
    }
    return priority;
}
std::tuple<int, Config_t, std::vector<JobID_t>> solve_greedy_ranges(const std::vector<JobInfo>& job_info, const std::vector<float>& priority, float coef = 1.f, float prio_coef=1.f) {
    std::vector<JobID_t> jobs(job_info.size());
    std::iota(jobs.begin(), jobs.end(), 0);

    double avg_duration = 0;
    for(auto i : job_info) avg_duration += i.duration;
    avg_duration /= job_info.size();

    std::sort(jobs.begin(), jobs.end(), [&](JobID_t a, JobID_t b) {
        return job_info[a].due - job_info[a].duration - (avg_duration * priority[a] * prio_coef) <
            job_info[b].due - job_info[b].duration - (avg_duration * priority[b] * prio_coef);
    });



    std::array<metric_t, NUM_MACHINES> machine_times;
    Config_t result;
    std::fill(machine_times.begin(), machine_times.end(), 0);

    std::vector<JobID_t> to_readd;
    int U = 0;
    for(auto job_id : jobs) {
        if(job_info[job_id].duration > avg_duration * coef) {
            to_readd.push_back(job_id);
            continue;
        }
        auto& job = job_info[job_id];
        size_t best_machine = -1;
        metric_t best_start = -1;
        for (int i = 0; i < NUM_MACHINES; i++) {
            metric_t start = std::max(job.ready, machine_times[i]);
            if(best_start == -1 || start < best_start) {
                best_start = start;
                best_machine = i;
            }
        }
        assert(best_start != -1);
        metric_t start = best_start;
        metric_t finish = start + job.duration;
        if(finish > job.due) {
            U += 1;
            to_readd.push_back(job_id);
        }else {
            machine_times[best_machine] = finish;
            result[best_machine].push_back(job_id);
        }
    }
    return std::make_tuple(U, result, to_readd);
}
//returns configuartion and the index at which the late jobs start
std::tuple<Config_t, int> solve(
    const std::vector<JobInfo>& job_info,
    const Config_t& greedy,
    std::vector<JobID_t> readd
) {
    std::vector<JobID_t> jobs(job_info.size());
    std::iota(jobs.begin(), jobs.end(), 0);
    IntervalTree<JobID_t> tree;
    for(auto j : jobs) {
        auto& job = job_info[j];
        tree.insert(j, job.ready, job.due - job.duration);
    }
    // most amount of ranges 
    auto ranges = search_ranges(job_info, jobs);
    
    Config_t result;
    int U;
    return std::make_tuple(result, U);
}
int findTardy(const std::vector<JobID_t>& machine, const std::vector<JobInfo>& job_info) {
    metric_t time = 0;
    int next_idx = 0;
    int U = 0;
    for(auto j : machine) {
        auto& job = job_info[j];
        time = std::max(time, job.ready);
        time += job.duration;
        if(time > job.due) {
            U += 1;
        }
    }
    return U;
}
//returns -1 if no free spots or idx for insert
int tryAddOne(std::vector<JobID_t> machine, JobID_t to_readd, const std::vector<JobInfo>& job_info) {
    metric_t time = 0;
    auto next_time = time;
    int next_idx = 0;
    auto & job = job_info[to_readd];
    while(next_time < job.ready && next_idx < machine.size()) {
        time = next_time;
        next_time = std::max(next_time, job_info[machine[next_idx]].ready);
        next_time += job_info[machine[next_idx]].duration;
        next_idx += 1;
    }
    while(next_idx < machine.size() && next_time < job.due) {
        machine.insert(machine.begin() + next_idx, to_readd);
        auto tardy = findTardy(machine, job_info);
        if(tardy == 0) {
            return next_idx;
        }
        machine.erase(machine.begin() + next_idx);
        if(next_idx == machine.size()) {
            break;
        }
        next_time = std::max(next_time, job_info[machine[next_idx]].ready);
        next_time += job_info[machine[next_idx]].duration;
        next_idx += 1;
    }
    return -1;
}
std::pair<int, JobID_t> tryAddOneWithPopping(std::vector<JobID_t> machine, JobID_t to_readd, const std::vector<JobInfo>& job_info) {
    metric_t time = 0;
    auto next_time = time;
    int next_idx = 0;
    auto & job = job_info[to_readd];
    while(next_time < job.ready && next_idx < machine.size()) {
        time = next_time;
        next_time = std::max(next_time, job_info[machine[next_idx]].ready);
        next_time += job_info[machine[next_idx]].duration;
        next_idx += 1;
    }
    while(next_idx < machine.size() && next_time < job.due) {
        machine.insert(machine.begin() + next_idx, to_readd);
        for(int to_pop = next_idx+1;next_idx < machine.size(); next_idx++) {
            auto elem = machine[to_pop];
            if(job_info[elem].duration <= job.duration)
                continue;

            machine.erase(machine.begin() + to_pop);
            auto tardy = findTardy(machine, job_info);
            if(tardy == 0) {
                return {next_idx, elem};
            }
            machine.insert(machine.begin() + to_pop, elem);
        }
        machine.erase(machine.begin() + next_idx);
        if(next_idx >= machine.size()) {
            break;
        }
        next_time = std::max(next_time, job_info[machine[next_idx]].ready);
        next_time += job_info[machine[next_idx]].duration;
        next_idx += 1;
    }
    return {-1, 0};
}

int calcU(const Config_t& sol, const std::vector<JobInfo>& job_info) {
    int U = 0;
    for(auto& m : sol) {
        U += findTardy(m, job_info);
    }
    return U;
}

Config_t tryReadd(Config_t& partial, std::vector<JobID_t> to_readd, const std::vector<JobInfo>& job_info) {
    auto prev_readd = to_readd;
    do {
        prev_readd = to_readd;
        std::set<JobID_t> added;
        for(auto& machine : partial) {
            for(auto j : to_readd) {
                if(added.find(j) != added.end())
                    continue;
                auto idx = tryAddOne(machine, j, job_info);
                if(idx != -1) {
                    machine.insert(machine.begin() + idx, j);
                    added.insert(j);
                    assert(findTardy(machine, job_info) == 0);
                }
            }
        }
        for(auto j : added) {
            auto itr = std::find(to_readd.begin(), to_readd.end(), j);
            to_readd.erase(itr);
        }
        added.clear();
        std::vector<JobID_t> to_readd_pp;
        for(auto& machine : partial) {
            for(auto j : to_readd) {
                if(added.find(j) != added.end())
                    continue;
                auto idx = tryAddOneWithPopping(machine, j, job_info);
                if(idx.first != -1) {
                    auto itr = std::find(machine.begin(), machine.end(), idx.second);
                    machine.erase(itr);
                    machine.insert(machine.begin() + idx.first, j);
                    added.insert(j);
                    to_readd_pp.push_back(idx.second);
                    assert(findTardy(machine, job_info) == 0);
                }
            }
        }
        for(auto j : added) {
            auto itr = std::find(to_readd.begin(), to_readd.end(), j);
            to_readd.erase(itr);
        }
        for(auto j : to_readd_pp) to_readd.push_back(j);
    }while(prev_readd != to_readd);
    for(auto j : to_readd) {
        partial[0].push_back(j);
    }
    return partial;
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

    if (signal(SIGALRM, alarm_handler) == SIG_ERR) {
        return 1;
    }
    alarm(sec_limit);

    if (!fin) {
        std::cerr << "Error opening file: " << argv[1] << std::endl;
        return 1;
    }
    
    std::vector<JobInfo> job_info;
    job_info = load_in_file(fin);
    int U_greedy = 0xffffff;
    Config_t greedy;
    Config_t full_greedy;
    std::vector<JobID_t> greedy_readd;

    auto prio = getPrio(job_info);
    for(float pcoef = -0.1; pcoef < 0.1; pcoef += 0.01) {
        for(float coef = 0.1; coef < 2; coef += 0.02) {
            auto [curU, curSol, readd] = solve_greedy_ranges(job_info, prio, coef, pcoef);
            // auto [curU, curSol, readd] = solve_greedy(job_info, coef);
            curU += readd.size();
            if(curU < U_greedy) {
                U_greedy = curU;
                greedy = curSol;
                greedy_readd = readd;
                for(auto j : readd) curSol[0].push_back(j);
                full_greedy = curSol;
            }
        }
    }
    for(auto m : greedy) {
        assert(findTardy(m, job_info) == 0);
    }
    auto solution = tryReadd(greedy, greedy_readd, job_info);
    auto U = calcU(solution, job_info);
    
    fout << U << "\n";
    for(auto m : solution) {
        for(auto job : m) {
            fout << job+1 << " ";
        }
        fout << "\n";
    }

    return 0;
}

