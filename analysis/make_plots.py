import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
import itertools
import math

ttsketch = ['seg', 'cs', 'srht']
vanilla = 'vanilla'

N = [
	'512', 
	# '1024', 
	# '2048',
]
# D = [4, 8, 16, 32, 64]
D = [4, 8]
V = range(20)
P = [2]
TRIALS = range(20)

K = [4**2, 8**2, 16**2, 32**2, 64**2, 128**2, 256**2]

for n, d, p in itertools.product(N, D, P):
	sketch_log_folder = '../tensorsketch/log/'
	vanilla_log_folder = '../../../data/vanillaLR/log/vanilla/'
	x_axis = list(K)
	for i in range(len(K)):
		x_axis[i] = int(math.sqrt(x_axis[i]))
	gamma = []
	sigma = []
	sigma_approx = []
	sigma_sketch = []
	for k in K:
		opt_time = 0
		opt_eval = 0
		for v in V:
			file_name = 'n_{}_d_{}_p_{}_v_{}.txt'.format(n, d, p, v)
			file_path = vanilla_log_folder + file_name
			file = open(file_path, 'r').read().split()
			opt_time += float(file[0])
			opt_eval += float(file[1])
		opt_time /= len(V)
		opt_eval /= len(V)

		sketch_time = 0
		approx_time = 0
		approx_eval = 0
		count = 0
		for v, trial in itertools.product(V, TRIALS):
			file_name = 'n_{}_d_{}_p_{}_v_{}_k_{}_trial_{}.txt'.format(n, d, p, v, k, trial)
			file_path = sketch_log_folder + file_name
			file = open(file_path, 'r').read().split()
			sketch_time += float(file[0])
			approx_time += float(file[1])
			approx_eval += float(file[2])
		sketch_time /= (len(V) * len(TRIALS))
		approx_time /= (len(V) * len(TRIALS))
		approx_eval /= (len(V) * len(TRIALS))

		gamma.append(approx_eval/opt_eval - 1)
		sigma.append((approx_time+sketch_time)/opt_time)
		sigma_approx.append(approx_time/opt_time)
		sigma_sketch.append(sketch_time/opt_time)
	
	gamma_ttsketch = {}
	sigma_ttsketch = {}
	sigma_approx_ttsketch = {}
	sigma_sketch_ttsketch = {}
	for prog in ttsketch:
		gamma_ttsketch[prog] = []
		sigma_ttsketch[prog] = []
		sigma_approx_ttsketch[prog] = []
		sigma_sketch_ttsketch[prog] = []
		sketch_log_folder = '../../../data/ttsketch/log/{}/'.format(prog)
		for k in x_axis:
			opt_time = 0
			opt_eval = 0
			for v in V:
				file_name = 'n_{}_d_{}_p_{}_v_{}.txt'.format(n, d, p, v)
				file_path = vanilla_log_folder + file_name
				file = open(file_path, 'r').read().split()
				opt_time += float(file[0])
				opt_eval += float(file[1])
			opt_time /= len(V)
			opt_eval /= len(V)

			sketch_time = 0
			approx_time = 0
			approx_eval = 0
			count = 0
			for v, trial in itertools.product(V, TRIALS):
				file_name = 'n_{}_d_{}_p_{}_v_{}_k_{}_trial_{}.txt'.format(n, d, p, v, k, trial)
				file_path = sketch_log_folder + file_name
				file = open(file_path, 'r').read().split()
				sketch_time += float(file[0])
				approx_time += float(file[1])
				approx_eval += float(file[2])
			sketch_time /= (len(V) * len(TRIALS))
			approx_time /= (len(V) * len(TRIALS))
			approx_eval /= (len(V) * len(TRIALS))

			gamma_ttsketch[prog].append(approx_eval/opt_eval - 1)
			sigma_ttsketch[prog].append((approx_time+sketch_time)/opt_time)
			sigma_approx_ttsketch[prog].append(approx_time/opt_time)
			sigma_sketch_ttsketch[prog].append(sketch_time/opt_time)



	plt.plot(x_axis, sigma, color='green', linestyle='dashed', linewidth = 1, 
	         marker='o', markerfacecolor='black', markersize=4, label='tensorsketch') 
	plt.plot(x_axis, sigma_ttsketch['cs'], color='red', linestyle='dashed', linewidth = 1, 
	         marker='o', markerfacecolor='black', markersize=4, label='cs') 
	plt.plot(x_axis, sigma_ttsketch['srht'], color='blue', linestyle='dashed', linewidth = 1, 
	         marker='o', markerfacecolor='black', markersize=4, label='srht') 
	plt.plot(x_axis, sigma_ttsketch['seg'], color='pink', linestyle='dashed', linewidth = 1, 
	         marker='o', markerfacecolor='black', markersize=4, label='seg') 
	plt.ylim(0, 3.0)
	plt.xlabel('K: sketch rows') 
	plt.ylabel('sigma: time/time') 
	plt.title('n_{}_d_{}_p_{}'.format(n, d, p)) 
	plt.legend()
	plt.savefig('./plots/' + 'sigma:n_{}_d_{}_p_{}.png'.format(n, d, p))
	plt.clf()

	plt.plot(x_axis, sigma_approx, color='green', linestyle='dashed', linewidth = 1, 
	         marker='o', markerfacecolor='black', markersize=4, label='tensorsketch') 
	plt.plot(x_axis, sigma_approx_ttsketch['cs'], color='red', linestyle='dashed', linewidth = 1, 
	         marker='o', markerfacecolor='black', markersize=4, label='cs') 
	plt.plot(x_axis, sigma_approx_ttsketch['srht'], color='blue', linestyle='dashed', linewidth = 1, 
	         marker='o', markerfacecolor='black', markersize=4, label='srht') 
	plt.plot(x_axis, sigma_approx_ttsketch['seg'], color='pink', linestyle='dashed', linewidth = 1, 
	         marker='o', markerfacecolor='black', markersize=4, label='seg') 
	plt.ylim(0, 3.0)
	plt.xlabel('K: sketch rows') 
	plt.ylabel('sigma_approx: time/time') 
	plt.title('n_{}_d_{}_p_{}'.format(n, d, p)) 
	plt.legend()
	plt.savefig('./plots/' + 'sigma_approx:n_{}_d_{}_p_{}.png'.format(n, d, p))
	plt.clf()

	plt.plot(x_axis, sigma_sketch, color='green', linestyle='dashed', linewidth = 1, 
	         marker='o', markerfacecolor='black', markersize=4, label='tensorsketch') 
	plt.plot(x_axis, sigma_sketch_ttsketch['cs'], color='red', linestyle='dashed', linewidth = 1, 
	         marker='o', markerfacecolor='black', markersize=4, label='cs') 
	plt.plot(x_axis, sigma_sketch_ttsketch['srht'], color='blue', linestyle='dashed', linewidth = 1, 
	         marker='o', markerfacecolor='black', markersize=4, label='srht') 
	plt.plot(x_axis, sigma_sketch_ttsketch['seg'], color='pink', linestyle='dashed', linewidth = 1, 
	         marker='o', markerfacecolor='black', markersize=4, label='seg') 
	plt.ylim(0, 3.0)
	plt.xlabel('K: sketch rows') 
	plt.ylabel('sigma_sketch: time/time') 
	plt.title('n_{}_d_{}_p_{}'.format(n, d, p)) 
	plt.legend()
	plt.savefig('./plots/' + 'sigma_sketch:n_{}_d_{}_p_{}.png'.format(n, d, p))
	plt.clf()

	# plt.plot(x_axis, gamma, color='green', linestyle='dashed', linewidth = 1, 
	#          marker='o', markerfacecolor='black', markersize=4, label='tensorsketch') 
	plt.plot(x_axis, gamma_ttsketch['cs'], color='red', linestyle='dashed', linewidth = 1, 
	         marker='o', markerfacecolor='black', markersize=4, label='cs') 
	# plt.plot(x_axis, gamma_ttsketch['srht'], color='blue', linestyle='dashed', linewidth = 1, 
	#          marker='o', markerfacecolor='black', markersize=4, label='srht') 
	# plt.plot(x_axis, gamma_ttsketch['seg'], color='pink', linestyle='dashed', linewidth = 1, 
	#          marker='o', markerfacecolor='black', markersize=4, label='seg') 
	plt.ylim(0, 0.2)
	plt.xlabel('K: sketch rows') 
	plt.ylabel('gamma: eval/eval - 1') 
	plt.title('n_{}_d_{}_p_{}'.format(n, d, p)) 
	plt.legend()
	plt.savefig('./plots/' + 'gamma:n_{}_d_{}_p_{}.png'.format(n, d, p))
	plt.clf()