import os
import sys
headers_path = os.environ['headers_path']
eigen_path = os.environ['eigen_path']
repo_home = os.environ['repo_home']

# arguments = sys.argv[1:]
# if len(arguments) == 0:
# 	print('*** ERROR: no commands ***')
# 	exit()

link_repo_home = '-I{}'.format(repo_home)
link_eigen = '-I{}'.format(eigen_path)

command = 'g++ -O0 main.cpp {} {} -std=c++11 -o .main.out'.format(link_eigen, link_repo_home)
os.system(command)
