#include <cassert>
#include <array>
#include <numeric>
#include <random>
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
    std::cout << "⏰ Alarm triggered!\n";
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
struct GAParams {
    int population_size = 500;
    int generations = 1000;
    float crossover_prob = 0.1;
    float mutation_prob = 0.05;
    int tournament_size = 3;
};
using Chromosome = std::vector<JobID_t>;
using Population = std::vector<Chromosome>;
int evaluate(const Chromosome &chrom, const std::vector<JobInfo>& job_info) {
    std::array<metric_t, NUM_MACHINES> machine_time;
    std::fill(machine_time.begin(), machine_time.end(), 0);
    int tardy_count = 0;

    for (int idx : chrom) {
        const JobInfo &job = job_info[idx];

        // Assign to earliest available machine respecting release time
        int best_machine = 0;
        long long start_time = std::max(machine_time[0], job.ready);
        for (int m = 1; m < NUM_MACHINES; ++m) {
            long long candidate = std::max(machine_time[m], job.ready);
            if (candidate < start_time) {
                start_time = candidate;
                best_machine = m;
            }
        }

        long long finish_time = start_time + job.duration;
        if (finish_time > job.due) tardy_count++;
        machine_time[best_machine] = finish_time;
    }

    return tardy_count;
}
Population initialize(int population_size, int num_jobs) {
    Population pop;
    pop.reserve(population_size);
    Chromosome base(num_jobs);
    std::iota(base.begin(), base.end(), 0);

    for (int i = 0; i < population_size; ++i) {
        std::shuffle(base.begin(), base.end(), rng);
        pop.push_back(base);
    }
    return pop;
}
Chromosome tournament_select(const Population &pop, const std::vector<int>& fitness, int k) {
    std::uniform_int_distribution<int> dist(0, pop.size() - 1);
    int best = dist(rng);
    for (int i = 1; i < k; ++i) {
        int cand = dist(rng);
        if (fitness[cand] < fitness[best]) best = cand;
    }
    return pop[best];
}
std::pair<Chromosome, Chromosome> crossover(const Chromosome &parent1, const Chromosome &parent2) {
    int n = parent1.size();
    std::uniform_int_distribution<int> dist(0, n-1);
    int l = dist(rng), r = dist(rng);
    if (l > r) std::swap(l, r);

    Chromosome child1(n, -1), child2(n, -1);

    // Copy segment from parent1
    for (int i = l; i <= r; ++i) child1[i] = parent1[i];
    for (int i = l; i <= r; ++i) child2[i] = parent2[i];

    auto fill_child = [&](const Chromosome &source, Chromosome &child){
        int pos = (r+1)%n;
        for (int i = 0; i < n; ++i) {
            int val = source[(r+1+i)%n];
            if (std::find(child.begin(), child.end(), val) == child.end()) {
                child[pos] = val;
                pos = (pos+1)%n;
            }
        }
    };

    fill_child(parent2, child1);
    fill_child(parent1, child2);

    return {child1, child2};
}
void mutate(Chromosome &chrom, float mutation_prob) {
    std::uniform_real_distribution<float> dist(0,1);
    if (dist(rng) < mutation_prob) {
        std::uniform_int_distribution<int> idx(0, chrom.size()-1);
        int a = idx(rng), b = idx(rng);
        std::swap(chrom[a], chrom[b]);
    }
}
Chromosome genetic_algorithm(const std::vector<JobInfo> &jobs, const GAParams &params = GAParams()) {
    int n = jobs.size();
    Population pop = initialize(params.population_size, n);

    std::vector<int> fitness(params.population_size);
    for (int i = 0; i < params.population_size; ++i)
        fitness[i] = evaluate(pop[i], jobs);

    Chromosome best = pop[0];
    int best_score = fitness[0];

    for (int gen = 0; gen < params.generations; ++gen) {
        Population new_pop;
        while (new_pop.size() < params.population_size) {
            // selection
            Chromosome p1 = tournament_select(pop, fitness, params.tournament_size);
            Chromosome p2 = tournament_select(pop, fitness, params.tournament_size);

            // crossover
            std::uniform_real_distribution<float> dist(0,1);
            std::pair<Chromosome, Chromosome> children;
            if (dist(rng) < params.crossover_prob) {
                children = crossover(p1, p2);
            } else {
                children = {p1, p2};
            }

            // mutation
            mutate(children.first, params.mutation_prob);
            mutate(children.second, params.mutation_prob);

            new_pop.push_back(children.first);
            if (new_pop.size() < params.population_size) new_pop.push_back(children.second);
        }

        // evaluate new population
        for (int i = 0; i < params.population_size; ++i) {
            fitness[i] = evaluate(new_pop[i], jobs);
            if (fitness[i] < best_score) {
                best_score = fitness[i];
                best = new_pop[i];
            }
        }

        pop = std::move(new_pop);
    }

    std::cout << "Best tardy jobs: " << best_score << "\n";
    return best;
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
std::vector<Range> search_ranges(const std::vector<JobInfo>& job_info, const std::vector<JobID_t>& jobs, int max_iters=10000) {
    metric_t max_due = 0;
    for(auto j : jobs) {
        max_due = std::max(max_due, job_info[j].due - job_info[j].duration);
    }
    float t_init = 0.5f;
    float t_min = 0.0f;
    float t_max = 1.0f;
    std::uniform_real_distribution<float> dist(-0.05f, 0.05f); // step size
    std::uniform_real_distribution<float> prob(0.0f, 1.0f);

    float current_t = t_init;
    long long current_score = get_ranges(job_info, jobs, (max_due * current_t)).size();

    float temperature = 1.0f;
    float cooling_rate = 0.995f;

    float best_t = current_t;
    long long best_score = current_score;

    for (int iter = 0; iter < max_iters; ++iter) {
        // propose a new t by small random step
        float new_t = current_t + dist(rng);
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
    std::cerr << "best threshold: " << best_t * max_due << "\t" << best_score << "\n";
    return get_ranges(job_info, jobs, (max_due * best_t));
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

    GAParams params;
    params.crossover_prob *= 0.f;
    params.mutation_prob = 0.1f;
    params.tournament_size *= 10.f;
    genetic_algorithm(job_info, params);
    
    Config_t result;
    int U;
    return std::make_tuple(result, U);
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
        std::cerr << "Signal callback could not be set";
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
    for(float coef = 0.1; coef < 2; coef += 0.01) {
        auto [curU, curSol, readd] = solve_greedy(job_info, coef);
        curU += readd.size();
        if(curU < U_greedy) {
            U_greedy = curU;
            greedy = curSol;
            greedy_readd = readd;
            for(auto j : readd) curSol[0].push_back(j);
            full_greedy = curSol;
        }
    }

    auto [U_squeeze, solution_squeeze] = solve(job_info, greedy, greedy_readd);
    auto U = U_greedy;
    auto solution = full_greedy;
    
    fout << U << "\n";
    for(auto m : solution) {
        for(auto job : m) {
            fout << job+1 << " ";
        }
        fout << "\n";
    }

    return 0;
}

