import sys
import glob
import logging
import argparse
import matplotlib.pyplot as plt

def check_files(inp, out):
    first_line = next(inp).split()
    n = int(first_line[0])
    setuptime = int(first_line[1])
    batch_size = int(first_line[2])
    in_data = inp.readlines()
    job_info = {}
    fig, ax = plt.subplots(figsize=(64, 32), nrows=2)
    for i, line in enumerate(in_data):
        vals = [int(x) for x in line.split()]
        procesing_time = vals[0]
        ready_time = vals[-1]
        job_info[i+1] = (ready_time, procesing_time)
    jobs = list(job_info.items())
    jobs.sort(key=lambda x: x[1][0])
    for i in range(len(jobs)):
        j = jobs[i][1]
        ax[0].barh(i, j[1], left=j[0], height=0.4)


    cmax = int(next(out))
    num_batches = int(next(out))
    out_data = out.readlines()
    out_num_jobs = 0
    last_finish_time = 0

    for line in out_data:
        vals = [int(x) for x in line.split()]
        if len(vals) > batch_size:
            raise ValueError(f"Batch size of {batch_size} restricts the solution")
        start = 0
        duration = 0
        for i in vals:
            start = max(start, job_info[i][0])
            duration = max(duration, job_info[i][1])
        start = max(start, last_finish_time)
            # ax.text(start + duration/2, i, job, va='center', ha='center', color='white', fontsize=10)
        for i in vals:
            ax[1].barh(i, job_info[i][1], left=start, height=0.4)
        stop_time = start + duration + setuptime
        last_finish_time = stop_time
        out_num_jobs += len(vals)
    last_finish_time -= setuptime
    last_finish_time = max(last_finish_time, 0)


    ax[0].set_xlabel("Time")
    ax[0].sharex(ax[1])
    ax[0].set_title("Job distribution")
    ax[1].set_xlabel("Time")
    ax[1].set_title("Job Scheduling Gantt Chart")
    plt.grid(axis='x', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig("out.png")
with open(sys.argv[1], 'r') as inp, open(sys.argv[2], 'r') as out:
    check_files(inp, out)
