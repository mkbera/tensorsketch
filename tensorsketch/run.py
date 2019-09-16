import os
repo_home = os.environ['repo_home']


logfile = '.temp'
n = 5
d = 2
p = 2
k = 2
folder = '{}/data/n_{}_d_{}_p_{}/'.format(repo_home, n, d, p)

args = ''
args += '--folder {} --logfile {} '.format(folder, logfile)
args += '--n {} --d {} --p {} --k {} '.format(n, d, p, k)
command = './.main.out {} '.format(args)

os.system(command)