#!/usr/bin/env python

import sys
import glob
import logging
import argparse

logging.basicConfig()
logger = logging.getLogger('dev')
logger.setLevel(logging.DEBUG)

STOP_ON_ERROR = True
CHECK_CMAX = True
CHECK_K = True
CHECK_NEGATIVE = True
CHECK_REPEATS = True
CHECK_OOB = True

CHECK_BSIZE = True
CHECK_NUMJOBS = True
CHECK_PARAMS = True
def check_out_file(f):
    cmax = int(next(f))
    num_batches = int(next(f))
    data = f.readlines()
    if num_batches != len(data) and CHECK_K:
        raise ValueError("Invalid k - number of batches does not match")

    if CHECK_PARAMS and num_batches <= 0:
        raise ValueError(f"Invalid batch number: {num_batches}")
    csum = 0
    unique = set()

    max_id = 0
    num_jobs = 0
    for line in data:
        vals = [int(x) for x in line.split()]

        num_jobs += len(vals)
        for x in vals:
            if x in unique and CHECK_REPEATS:
                raise ValueError(f"Job {x} is repeated")
            unique.add(x)
            if x <= 0 and CHECK_OOB:
                raise ValueError(f"Job id must be greater than 0, but is {x}")

        max_v = max(vals)
        max_id = max(max_v, max_id)
    if max_id != num_jobs and CHECK_OOB:
        raise ValueError(f"Job {max_id} is too high of an id")
def check_in_file(inp):
    first_line = next(inp).split()
    n = int(first_line[0])
    setuptime = int(first_line[1])
    batch_size = int(first_line[2])
    if CHECK_PARAMS and (n <= 0 or setuptime < 0 or batch_size <= 0):
        raise ValueError(f"Invalid input parameters")

    in_data = inp.readlines()
    job_info = {}
    for i, line in enumerate(in_data):
        vals = [int(x) for x in line.split()]
        procesing_time = vals[0]
        ready_time = vals[-1]
        if CHECK_PARAMS and procesing_time < 0:
            raise ValueError("Processing time cannot be less than zero")
        if CHECK_PARAMS and ready_time < 0:
            raise ValueError("Ready time cannot be less than zero")

def check_files(inp, out):
    first_line = next(inp).split()
    n = int(first_line[0])
    setuptime = int(first_line[1])
    batch_size = int(first_line[2])
    in_data = inp.readlines()
    job_info = {}
    for i, line in enumerate(in_data):
        vals = [int(x) for x in line.split()]
        procesing_time = vals[0]
        ready_time = vals[-1]
        job_info[i+1] = (ready_time, procesing_time)


    cmax = int(next(out))
    num_batches = int(next(out))
    out_data = out.readlines()
    out_num_jobs = 0
    last_finish_time = 0
    for line in out_data:
        vals = [int(x) for x in line.split()]
        if CHECK_BSIZE and len(vals) > batch_size:
            raise ValueError(f"Batch size of {batch_size} restricts the solution")
        start = 0
        for i in vals:
            start = max(start, job_info[i][0])
        start = max(start, last_finish_time)
        duration = 0
        for i in vals:
            duration = max(duration, job_info[i][1])
        stop_time = start + duration + setuptime
        last_finish_time = stop_time
        out_num_jobs += len(vals)
    last_finish_time -= setuptime
    last_finish_time = max(last_finish_time, 0)

    if CHECK_NUMJOBS and out_num_jobs != n:
        raise ValueError(f"Number of jobs does not match {out_num_jobs}!")
    if CHECK_CMAX and last_finish_time != cmax:
        logger.warning(f"Cmax_read: {cmax}, Cmax_calculated: {last_finish_time}")
        raise ValueError(f"calculated cmax does not match {last_finish_time} != {cmax}")
    return cmax


def main(in_file, out_file):
    failed_already = None
    try:
        logger.debug(f"Processing {in_file},{out_file}:")
        with open(out_file, 'r') as f:
            check_out_file(f)
        with open(in_file, 'r') as f:
            check_in_file(f)
        with open(in_file, 'r') as inp, open(out_file, 'r') as out:
            cmax = check_files(inp, out)
            print(cmax)
    except ValueError as e:
        if STOP_ON_ERROR:
            raise ValueError(f"in {out_file}: {e}")
        else:
            if not failed_already:
                failed_already = ""
            failed_already += f"in {out_file}: {e}\n"


    if failed_already:
        raise ValueError(failed_already)
    return 0

parser = argparse.ArgumentParser(
            prog='1B-checker',
            description='This program generates scenarions for the 1B problem.')

parser.add_argument('in_file')
parser.add_argument('out_file')
parser.add_argument('-s', '--stop_on_error', action='store_true', required=False)
parser.add_argument('-c', '--check_cmax', action='store_true', required=False)
parser.add_argument('-k', '--check_k', action='store_true', required=False)
parser.add_argument('-n', '--check_neg', action='store_true', required=False)
parser.add_argument('-r', '--check_rep', action='store_true', required=False)
parser.add_argument('-b', '--check_bound', action='store_true', required=False)


if __name__ == "__main__":
    args = parser.parse_args()
    if args.stop_on_error:
        STOP_ON_ERROR = True 
    if args.check_cmax:
        CHECK_CMAX = True 
    if args.check_k:
        CHECK_K = True 
    if args.check_neg:
        CHECK_NEGATIVE = True 
    if args.check_rep:
        CHECK_REPEATS = True 
    if args.check_bound:
        CHECK_REPEATS = True 
    try:
        sys.exit(main(args.in_file, args.out_file))
    except Exception as e:
        logger.error(f"Err: {e}")
        sys.exit(1)
    sys.exit(0)
