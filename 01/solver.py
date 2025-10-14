import sys
import glob
import logging
import random
import time
import argparse
from collections import namedtuple
JobInfo = namedtuple('JobInfo', ['ready', 'duration'])

def finish_time(job):
    return job.ready + job.duration


def check_in_file(inp):
    first_line = next(inp).split()
    n = int(first_line[0])
    setuptime = int(first_line[1])
    batch_size = int(first_line[2])
    if (n <= 0 or setuptime < 0 or batch_size <= 0):
        raise ValueError(f"Invalid input parameters")
    print(f"Batch size: {batch_size}", file=sys.stderr)

    in_data = inp.readlines()
    job_info = [None] * n
    for i, line in enumerate(in_data):
        vals = [int(x) for x in line.split()]
        procesing_time = vals[0]
        ready_time = vals[-1]
        if procesing_time < 0:
            raise ValueError("Processing time cannot be less than zero")
        if ready_time < 0:
            raise ValueError("Ready time cannot be less than zero")
        job_info[i] = JobInfo(ready_time, procesing_time)
    return setuptime, batch_size, job_info

def solve_greedy(setuptime, batch_size, job_info):
    jobs = [x for x in range(len(job_info))]
    jobs.sort(key=lambda x: -job_info[x].ready)
    time = job_info[jobs[-batch_size]].ready
    print(f"Time as the start: {time}", file=sys.stderr)
    batches = []
    while len(jobs) != 0:
        if batches == [] or batches[-1] != []:
            batches.append([])
        left = len(jobs) - 1
        while left >= 0 and job_info[jobs[left]].ready <= time:
            left -= 1
        potential = jobs[left+1:]
        potential.sort(key=lambda x: job_info[x].duration)
        new_time = time
        for job in potential[:batch_size]:
            if job_info[job].ready > time:
                break
            new_time = max(time + job_info[job].duration, new_time)
            batches[-1].append(job)
            jobs.remove(job)
        time = new_time + setuptime

    if batches[-1] == []:
        batches.pop()
    return batches


def calc_cmax(setuptime, batch_size, batches, job_info):
    cmax = 0
    for batch in batches:
        ready = 0
        duration = 0
        for job in batch:
            ready = max(ready, job_info[job].ready)
            duration = max(duration, job_info[job].duration)
        cmax = max(ready, cmax) 
        cmax += duration
        cmax += setuptime
    cmax -= setuptime
    return cmax

def perform_swap(solution, batch_id, id1, id2):
    b1 = batch_id[id1]
    b2 = batch_id[id2]
    solution[b1].remove(id1)
    solution[b2].remove(id2)

    solution[b1].append(id2)
    solution[b2].append(id1)
    return solution

def good_swap(i1, i2, tabu_list, min_iter=-1):
    if i1 == i2:
        return False
    id1 = min(i1, i2)
    id2 = max(i1, i2)
    if tabu_list[id1] is not None and tabu_list[id1][id2] is not None:
        if tabu_list[id1][id2] > min_iter:
            return False
        else:
            return True
    return True

def add_to_tabu(tabu, id1, id2, n, iter=1):
    if tabu[id1] is None:
        tabu[id1] = [None] * n

    tabu[id1][id2] = iter


def beam_search(setuptime, batch_size, batches, job_info, num_beams=7, num_rand=2, time_limit=5):
    start = time.time()
    n = len(job_info)
    iter_size = n
    global_tdur = 100
    c = calc_cmax(setuptime, batch_size, batches, job_info)

    global_tabu_list = [None] * n
    winners = [batches]
    last_swaps = []
    global_iter = 0
    bast_cmax_so_far = float('inf')
    best_winner = None
    while True:
        if time.time() - start > time_limit:
            break
        global_iter += 1
        swaps = []
        for wid, solution in enumerate(winners):
            batch_id = [None] * n
            for i, b in enumerate(solution):
                for job in b:
                    batch_id[job] = i
            local_tabu_list = [None] * n
            for local_iter in range(iter_size):
                i1 = random.randrange(0, n)
                i2 = i1
                while batch_id[i1] == batch_id[i2] or not good_swap(i1, i2, local_tabu_list) or not good_swap(i1, i2, global_tabu_list, global_iter-global_tdur):
                    i2 = random.randrange(0, n)
                id1 = min(i1, i2)
                id2 = max(i1, i2)
                add_to_tabu(local_tabu_list, id1, id2, n)

                temp = perform_swap([r[:] for r in solution], batch_id, id1, id2)

                cmax = calc_cmax(setuptime, batch_size, temp, job_info)
                swaps.append((cmax, temp, id1, id2))
        # calc winners
        swaps = swaps + last_swaps

        swaps.sort(key=lambda i: i[0])
        last_swaps = []
        winners = []
        max_win_id = min(len(swaps), num_beams)
        win_idx = [x for x in range(max_win_id)]
        rand_idx = [random.randrange(max_win_id, len(swaps)) for _ in range(num_rand)]
        choosen = win_idx + rand_idx
        for i in choosen:
            last_swaps.append(swaps[i])
            winners.append(swaps[i][1])
            id1 = swaps[i][2]
            id2 = swaps[i][3]
            add_to_tabu(global_tabu_list, id1, id2, n, global_iter)
        if swaps[0][0] < bast_cmax_so_far:
            bast_cmax_so_far = swaps[0][0]
            best_winner = swaps[0][1]

        print(f"Best of iteration {global_iter}: {swaps[0][0]}", file=sys.stderr)
    return best_winner

random.seed(42)
with open(sys.argv[1], 'r') as f:
    setup, batch_s, job_info = check_in_file(f)
solution = solve_greedy(setup, batch_s, job_info)
solution = beam_search(setup, batch_s, solution, job_info)

print(f"0") #cmax
print(len(solution)) #num-batches
for b in solution: #batches
    for j in b: #solution
        print(j+1, end=" ")
    print("")
print(calc_cmax(setup, batch_s, solution, job_info), file=sys.stderr)

