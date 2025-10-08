#!/usr/bin/env python

import argparse
import sys
import glob
import logging
import math
import random

_box_muller_spare = None

def normal_rand(mu=0.0, sigma=1.0):
    global _box_muller_spare
    if sigma < 0:
        raise ValueError("sigma must be non-negative")
    if _box_muller_spare is not None:
        z = _box_muller_spare
        _box_muller_spare = None
        return mu + sigma * z
    u1 = 0.0
    while u1 == 0.0:
        u1 = random.random()
    u2 = random.random()
    r = math.sqrt(-2.0 * math.log(u1))
    theta = 2.0 * math.pi * u2
    z0 = r * math.cos(theta)
    z1 = r * math.sin(theta)
    _box_muller_spare = z1
    return mu + sigma * z0
def gamma_rand(alpha, theta=1.0):
    if alpha <= 0 or theta <= 0:
        raise ValueError("alpha and theta must be positive")
    if alpha < 1.0:
        while True:
            u = random.random()
            b = (math.e + alpha) / math.e
            p = b * u
            if p <= 1.0:
                x = p ** (1.0 / alpha)
            else:
                x = -math.log((b - p) / alpha)
            u2 = random.random()
            if p <= 1.0:
                if u2 <= math.exp(-x):
                    return theta * x
            else:
                if u2 <= x ** (alpha - 1):
                    return theta * x
    else:
        d = alpha - 1.0 / 3.0
        c = 1.0 / math.sqrt(9.0 * d)
        while True:
            while True:
                u1 = random.random()
                u2 = random.random()
                v = math.sqrt(-2.0 * math.log(u1)) * math.cos(2.0 * math.pi * u2)
                x = 1.0 + c * v
                if x > 0:
                    break
            x3 = x * x * x
            u = random.random()
            if u < 1.0 - 0.0331 * (v ** 4):
                return theta * d * x3
            if math.log(u) < 0.5 * v * v + d * (1 - x3 + math.log(x3)):
                return theta * d * x3

def exp_rand(lmbda=1.0):
    if lmbda <= 0:
        raise ValueError("lambda must be positive")
    u = 0.0
    while u == 0.0:
        u = random.random()
    return -math.log(u) / lmbda

def uniform_rand(lmbda=1.0):
    if lmbda <= 0:
        raise ValueError("lambda must be positive")
    return random.random() * lmbda

def beta_rand(alpha, beta):
    if alpha <= 0 or beta <= 0:
        raise ValueError("alpha and beta must be positive")

    x = gamma_rand(alpha, 1.0)
    y = gamma_rand(beta, 1.0)
    return x / (x + y)
def normalize(data):
    if not data:
        return []
    dmin = min(data)
    dmax = max(data)
    if dmax == dmin:
        return [0.0 for _ in data]
    return [(x - dmin) / (dmax - dmin) for x in data]

mu, sigma = 0, 1
alpha, theta = 2.0, 2.0
alpha_b, beta_b = 2.0, 5.0
lmbda = 1.5

def gen_normal(N):
    return [normal_rand(mu, sigma) for _ in range(N)]
def gen_gamma(N):
    return [gamma_rand(alpha, theta) for _ in range(N)]
def gen_exp(N):
    return [exp_rand(lmbda) for _ in range(N)]
def gen_beta(N):
    return [beta_rand(alpha_b, beta_b) for _ in range(N)]
def gen_uniform(N):
    return [random.random() for _ in range(N)]

gen_func = [gen_normal, gen_gamma, gen_exp, gen_beta, gen_uniform]

BIG_NUMBER = 100
def generate(N, f):
    gen = gen_func[random.randrange(len(gen_func))]
    ready_times = sorted(normalize(gen(N)))
    ready_times = [x * BIG_NUMBER for x in ready_times]
    avg_gap = []
    for i in range(N - 1):
        avg_gap.append(ready_times[i + 1] - ready_times[i])
    median_gap = sorted(avg_gap)[len(avg_gap)//2]
    avg = sum(ready_times) / N

    choosen_batch_size = random.randrange(N//10, N//5)
    choosen_cooloff = int(choosen_batch_size * 0.1 * math.sqrt(avg))
    f.write(f"{N} {choosen_cooloff} {choosen_batch_size}\n")

    durations = []
    duration_mean = median_gap * (choosen_batch_size)
    duration_sd = median_gap * (choosen_batch_size-1) * 0.25
    for _ in range(N):
        durations.append( normal_rand(duration_mean, duration_sd) )

    for dur, ready in zip(durations, ready_times):
        f.write(f"{int(dur+1)} {int(ready+1)}\n")

parser = argparse.ArgumentParser(
            prog='1B-generator',
            description='This program generates scenarions for the 1B problem.')

parser.add_argument('-n', '--num_jobs', type=int, required=True, help='number of jobs that will be generated.')      # option that takes a value
parser.add_argument('-m', '--max_num', type=int, help='upper bound of the jobs ready time.')
parser.add_argument('-o', '--file', required=True, help='output file that will store generated scenario.')
args = parser.parse_args()

if args.max_num:
    BIG_NUMBER = args.max_num
else:
    BIG_NUMBER = args.num_jobs * 3
with open(args.file, "w") as f:
    generate(args.num_jobs, f)


