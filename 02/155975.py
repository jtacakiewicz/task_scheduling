#!/usr/bin/env python

import sys
import glob
import logging
import argparse

logging.basicConfig()
logger = logging.getLogger('dev')
logger.setLevel(logging.DEBUG)

MACHINE_COUNT = 5

def validate_in(inp):
    first_line = next(inp).split()
    n = int(first_line[0])
    jobs = {}
    for i, line in enumerate(inp):
        params=line.split()
        dur = int(params[0])
        rdy = int(params[1])
        due = int(params[2])
        if dur <= 0:
            logger.warning(f"Job {i+1} has invalid duration: {dur}")
        if rdy < 0:
            logger.warning(f"Job {i+1} has invalid ready time: {rdy}")
        if rdy+dur > due:
            logger.warning(f"Job {i+1} is undoable in these requirements")
        jobs[i+1] = (dur, rdy, due)

    if len(jobs) != n:
        logger.error(f"Input file did not describe all jobs ({len(jobs)})")

    return n, jobs

def solve_greedy(n, job_info, coef = 1):
    jobs = []
    for id, info in job_info.items():
        jobs.append((id, *info))
    avg_dur = 0
    durations = map(lambda x: x[1], jobs)
    avg_dur = sum(durations) / n

    jobs.sort(key=lambda x: x[-1])
    machine_time = [0] * MACHINE_COUNT
    machine_jobs = [[] for _ in range(MACHINE_COUNT)]
    U = 0
    completion_times = [0] * (n + 1)

    to_reapply = []
    for job_id, dur, rdy, due in jobs:
        if dur > avg_dur * coef:
            to_reapply.append(job_id)
            continue;
        best_machine = None
        best_start = None

        for m in range(MACHINE_COUNT):
            start = max(rdy, machine_time[m])
            if best_start is None or start < best_start:
                best_start = start
                best_machine = m
        
        start = best_start
        finish = start + dur

        machine_time[best_machine] = finish
        machine_jobs[best_machine].append(job_id)
        completion_times[job_id] = finish

        if finish > due:
            U += 1
    for id in to_reapply:
        machine_jobs[0].append(id)
        U += 1
    return U, machine_jobs


parser = argparse.ArgumentParser(
            prog='2B-solver',
            description='This program solves greedly the 2B problem.')

parser.add_argument('in_file', help='input file.')
parser.add_argument('out_file', help='output file.')

args = parser.parse_args()
with open(args.in_file, "r") as inp:
    n, job_info = validate_in(inp)
    best_sol = None
    best_U = None
    for coef_int in range(0, 5000):
        coef = coef_int / 1000
        U, sol = solve_greedy(n, job_info, coef)
        if best_U is None or U < best_U:
            best_U = U
            best_sol = sol

    with open(args.out_file, "w") as out:
        out.write(str(best_U) + '\n')
        for m in range(MACHINE_COUNT):
            out.write(" ".join(map(str, best_sol[m])) + '\n')
