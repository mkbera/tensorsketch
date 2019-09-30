import os
import sys
headers_path = os.environ['headers_path']
eigen_path = os.environ['eigen_path']
repo_home = os.environ['repo_home']

arguments = sys.argv[1:]
if len(arguments) == 0:
	print('*** ERROR: no commands ***')
	exit()
arg = arguments[0]
link_repo_home = '-I{}'.format(repo_home)
link_eigen = '-I{}'.format(eigen_path)

if arg == '--ts':
	command = 'g++ -O3 ts.cpp {} {} -std=c++11 -o .ts.out'.format(link_eigen, link_repo_home)
if arg == '--ts2':
	command = 'g++ -O3 ts2.cpp {} {} -std=c++11 -o .ts2.out'.format(link_eigen, link_repo_home)
os.system(command)
