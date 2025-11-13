#!/usr/bin/env python

import argparse
import sys
import glob
import logging
import math
import random
import numpy as np

JOB_DURATION = 40
NUM_MACHINES=5

def generate_cheat(N, f, ratio=0.5):
    avg_duration_og = JOB_DURATION
    avg_spare_og = avg_duration_og*0.75
    std_og = 0.3

    avg_duration_surplus = JOB_DURATION*1.2
    avg_spare_surplus = avg_duration_surplus*0.6
    std_surplus = 0.5
    machN = N // NUM_MACHINES

    surplus = int(machN * ratio)
    og = machN - surplus
    # processtime, readytime, duedate
    jobs = []
    for j in range(NUM_MACHINES):
        last_end = 0
        min_ready = 0
        max_ready = 0
        for _ in range(og):
            duration    = int(max(np.random.normal(avg_duration_og, avg_duration_og*std_og), 1))
            ready_spare = int(np.random.normal(avg_spare_og, avg_spare_og*std_og*0.1))
            due_spare   = int(np.random.normal(avg_spare_og, avg_spare_og*std_og*0.1))
            jobs.append((duration, max(last_end-ready_spare, 0), last_end+duration+due_spare))
            last_end = last_end + duration
            max_ready = last_end
        print(max_ready)

        for _ in range(surplus):
            readytime = int(np.random.uniform(min_ready, max_ready-avg_spare_surplus))
            duration    = int(max(np.random.normal(avg_duration_surplus, avg_duration_surplus*std_surplus), 1))
            ready_spare = int(np.random.normal(avg_spare_surplus, avg_spare_surplus*std_surplus*0.1))
            due_spare   = int(np.random.normal(avg_spare_surplus, avg_spare_surplus*std_surplus*0.1))
            jobs.append((duration, max(readytime-ready_spare, 0), readytime+duration+due_spare))
    np.random.shuffle(jobs)

    f.write(f"{N}\n")
    for dur, rdy, due in jobs:
        f.write(f"{dur} {rdy} {due}\n")

parser = argparse.ArgumentParser(
            prog='2B-generator',
            description='This program generates scenarions for the 2B problem.')

parser.add_argument('-n', '--num_jobs', type=int, required=True, help='number of jobs that will be generated.')      # option that takes a value
parser.add_argument('-o', '--file', required=True, help='output file that will store generated scenario.')
parser.add_argument('-j', '--job_duration', type=int, help='avg job duration.')
parser.add_argument('-p', '--num_machines', type=int, help='number of machines.')

args = parser.parse_args()

if args.num_machines:
    NUM_MACHINES = args.num_machines

if args.job_duration:
    JOB_DURATION = args.job_duration

with open(args.file, "w") as f:
    generate_cheat(args.num_jobs, f)
