import os
import os
import itertools
import sys
import random
data_path = os.environ['data_path']
repo_home = os.environ['repo_home']
utilities_path = os.environ['utilities_path']
sys.path.append(utilities_path)
from parallel_prog import run_parallel_progs_in_threads as RPPIT


repo_home = os.environ['repo_home']

# N = [512]
# # D = [4, 8, 16, 32, 64]
# D = [16, 32, 64]
# V = range(20)
# P = [2]
# K = [256**2, 128**2, 64**2, 32**2, 16**2, 8**2, 4**2]
# # K = [256**2, 128**2,]
# TRIALS = range(20)

N = [512]
D = [4]
V = [0]
P = [2]
K = [ 256**2]
TRIALS = [0]

# prog = 'ts2'
prog = 'ts'
cmd_list = []
for n, d, p, k, v, trial in itertools.product(N, D, P, K, V, TRIALS):
    folder_path = data_path + 'n_{}_d_{}_p_{}_v_{}/'.format(n, d, p, v)
    # log_file = repo_home + '.temp/log.txt'
    log_file = './log/' + 'n_{}_d_{}_p_{}_v_{}_k_{}_trial_{}.txt'.format(n, d, p, v, k, trial)
    try:
    	os.makedirs('./log/')
    except:
    	pass
    args = '--folder {} '.format(folder_path)
    # args += '--logfile {} '.format(log_file)
    args += '--logfile {} '.format('.temp')
    args += '--n {} --d {} --p {} --k {} '.format(n, d, p, k)
    seed = random.randint(1, 1000000)
    args += '--seed {} '.format(seed)
    command = './.{}.out {}'.format(prog, args)
    cmd_list.append(command)
RPPIT(cmd_list, 10)
