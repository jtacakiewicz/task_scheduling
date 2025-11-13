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

def validate_out(out, n):
    first_line = next(out).split()
    U = int(first_line[0])
    machines = []
    jobs_proc = set()
    for line in out:
        machines.append([])
        for job in line.split():
            if job in jobs_proc:
                logger.warning(f"Job {job} has been used twice")
            jobs_proc.add(job)
            machines[-1].append(job)
    if len(machines) > MACHINE_COUNT:
        logger.error(f"Too many machines: {len(machines)}")
    if len(jobs_proc) != n:
        logger.error(f"Incorrect amount of jobs used: {len(jobs_proc)}/{n}")
    return U, machines

def validate_solution(solution, job_info, U_exp):
    U = 0
    for machine in solution:
        last_end = 0
        for job in machine:
            job = int(job)
            if job not in job_info:
                logger.error(f"Job used in solution does not exist: {job}")
                continue;

            dur = job_info[job][0]
            rdy = job_info[job][1]
            due = job_info[job][2]

            # wait for ready time
            if rdy > last_end:
                last_end = rdy

            this_end = last_end + dur
            if this_end > due:
                U += 1
    if U != U_exp:
        logger.warning(f"Calculated U is different: {U} != {U_exp}")

    return U





parser = argparse.ArgumentParser(
            prog='2B-validator',
            description='This program validates solution of the 2B problem.')

parser.add_argument('in_file', help='input file.')
parser.add_argument('out_file', help='output file.')

args = parser.parse_args()
with open(args.in_file, "r") as inp:
    with open(args.out_file, "r") as out:
        n, job_info = validate_in(inp)
        U, sol = validate_out(out, n)
        validate_solution(sol, job_info, U)
