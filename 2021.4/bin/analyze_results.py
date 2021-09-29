""" 
This script analyzes the partition results 
"""
import networkx as nx
import os 
import matplotlib.pyplot as plt 
import matplotlib.ticker as mticker
import csv
import numpy as np
import genetic_partition_test as gp
import itertools
from networkx.drawing.nx_agraph import graphviz_layout
import matplotlib.image as mpimg
import seaborn as sns

def get_node_number (edgelist):
	G = nx.read_edgelist (edgelist, nodetype = str, create_using=nx.DiGraph())
	return len(list(G.nodes()))

def load_graph (edgelist):
	G = nx.read_edgelist (edgelist, nodetype = str, create_using=nx.DiGraph())
	return G



def load_data (filename):
	"""Load data file"""
	data = {}
	data_reader = csv.reader(open(filename, 'rU'), delimiter='\t')
	# Ignore header
	header = next(data_reader)
	# Process each line
	sol = 1
	for row in data_reader:
		if len(row) == len(header):
			sample = sol
			sample_data = {}
			for el_idx, el in enumerate(header):
				sample_data[el] = row[el_idx]
			data[sample] = sample_data
			sol += 1
	return data

def load_sol (filename):
	"""Load best_solns file"""
	data = {}
	data_reader = csv.reader(open(filename, 'rU'), delimiter='\t')
	# Ignore header
	header = next(data_reader)
	# Process each line
	for row in data_reader:
		if len(row) == len(header):
			sample = row[1]
			sample_data = {}
			for el_idx, el in enumerate(header[2:],2):
				sample_data[el] = row[el_idx]
			data[sample] = sample_data
	return data

def plot_network (G, outfig):

	pos = nx.kamada_kawai_layout(G)
	plt.figure(num=None, figsize=(8,8), dpi=80)
	nx.draw(
	    G,
	    pos=pos,
	    horizontalalignment='left',
	    verticalalignment='bottom',
	    node_color='powderblue'
	)
	plt.savefig(outfig, dpi=200)
	plt.show()


def count_nodes (indir1, indir2):
	""" count the number of nodes in graphs """
	count1, count2 = [], []
	for benchmark in os.listdir(indir1):
		if benchmark.startswith('.') == False:	
			edgelist1    = PATH + 'runs/benchmark/4-input-boolean-circuits/' + benchmark + '/DAG.edgelist'
			node_number1 = get_node_number (edgelist1)
			count1.append( node_number1 )
	for benchmark in os.listdir(indir2):
		if benchmark.startswith('.') == False:	
			edgelist2    = PATH + 'runs/benchmark/5-input-boolean-circuits/' + benchmark + '/DAG.edgelist'
			node_number2 = get_node_number (edgelist2)
			count2.append( node_number2 )

	fig = plt.figure(figsize=(5,5)) 
	ax  = fig.add_subplot(111)
	colors = ['orange', 'blue']
	labels = ['4-input Boolean circuits', '5-input Boolean circuits']
	plt.hist([count1, count2], 10, histtype='step', stacked=False, fill=False, label=labels)
	# ax.hist(count1, ec='k', histtype='step', fill=False, facecolor='wheat')
	# ax.hist(count2, ec='k', histtype='step', fill=False, facecolor='g')
	ax.set_xlabel('Number of vertices')
	ax.set_ylabel('Count')
	ax.set_xlim([0, 50])
	plt.legend(loc='upper right', frameon=False)
	plt.gca().yaxis.set_major_locator(mticker.MultipleLocator(4))
	plt.gca().xaxis.set_major_locator(mticker.MultipleLocator(5))
	for spine in ['top', 'right']:
		ax.spines[spine].set_linewidth(0)
	plt.savefig('Node distribution.pdf', dpi=200)
	plt.show()

def partition_stats ():
	ori_hc_valid_count, ori_lc_valid_count, hc_valid_count, lc_valid_count = 0, 0, 0, 0
	bm_path = PATH + 'runs/benchmark/5-input-boolean-circuits/'
	rs_path = PATH + 'runs/results/5-input-boolean-circuits/'

	for benchmark in os.listdir(bm_path):
		# determine whether original partition satisfy hc or lc
		if benchmark.startswith('.') == False:	
			# check whether original partition satisfy lc/hc 
			ori_hc_valid_bm, ori_lc_valid_bm = 0, 0
			edgelist = bm_path + benchmark + '/DAG.edgelist'
			G = load_graph (edgelist)
			in_nodes, out_nodes, nonprimitives  = gp.get_nonprimitive_nodes (G)
			G_primitive = gp.get_G_primitive (G, nonprimitives)
			for npart in os.listdir(rs_path + benchmark + '/nparts/'):
				if npart.isdigit():
					part_sol = rs_path + benchmark + '/nparts/' + npart + '/part_solns.txt'
					cut, partDict = gp.load_metis_part_sol (part_sol)
					part_opt = (cut, [gp.get_part(partDict, n) for n in G_primitive.nodes()])
					matrix, partG = gp.partition_matrix (G_primitive, part_opt)
					loop_free = gp.check_cycles(partG)
					motif_allowed_hc = gp.check_motif_allowed(matrix, '2,2')
					motif_allowed_lc = gp.check_motif_allowed(matrix, '4')
					if loop_free and motif_allowed_hc: 
						ori_hc_valid_bm += 1
					if loop_free and motif_allowed_lc: 
						ori_lc_valid_bm += 1

			# check whehter optimized solution satisfy lc/hc
			best_sol = rs_path + benchmark + '/nparts/best_solns.txt'

			if os.path.exists(best_sol):
				data = load_data (best_sol)
				hc_valid_bm, lc_valid_bm = 0, 0
				for sol in data:
					constraint = data[sol]['Constraint']
					if constraint == 'hc': 
						hc_valid_bm += 1
					else: 
						lc_valid_bm += 1
				
				if ori_hc_valid_bm > 0:
					ori_hc_valid_count += 1
				if ori_lc_valid_bm > 0:
					ori_lc_valid_count += 1
				if hc_valid_bm > 0:
					hc_valid_count += 1
				if lc_valid_bm > 0: 
					lc_valid_count += 1

			else: 
				print('best sol for benchmark', benchmark, 'does not exist')

	f_out = open (PATH+'runs/results/analysis/Constraint stats-5-input.txt', 'w')
	f_out.write('Original Partition\thc\t'+str(ori_hc_valid_count)+'\n')
	f_out.write('Original Partition\tlc\t'+str(ori_lc_valid_count)+'\n')
	f_out.write('Optimized Partition\thc\t'+str(hc_valid_count)+'\n')
	f_out.write('Optimized Partition\tlc\t'+str(lc_valid_count)+'\n')


def compile_best_solutions ():

	bm_path = PATH + 'runs/benchmark/5-input-boolean-circuits/'
	rs_path = PATH + 'runs/results/5-input-boolean-circuits/'

	f_out = open (PATH+'runs/results/5-input-boolean-circuits/best_solns.txt', 'w')
	f_out.write('\t'.join(['Category', 'Circuit', 'Npart', 'Sol', 'Nodes','Constraint','Valid Motif_METIS','Valid Motif_Optimized','Cycle Free_METIS','Cycle Free_Optimized','T_Metis','T_Optimized','cut_Metis','cut_Optimized'])+'\n')
	for benchmark in os.listdir(bm_path):
		if benchmark.startswith('.') == False:	
			best_sol = rs_path + benchmark + '/nparts/best_solns.txt'
			if os.path.exists(best_sol):
				data = load_data (best_sol)
				if data != {}: 
					best_sol = data[1]
					hc_sol = [sol for sol in data.keys() if data[sol]['Constraint'] == 'hc']
					lc_sol = [sol for sol in data.keys() if data[sol]['Constraint'] == 'lc']
					if hc_sol != []:
						# first choose from hc sols
						minT = min([data[sol]['T_Optimized'] for sol in hc_sol])
						minN = min([data[sol]['Npart'] for sol in hc_sol if data[sol]['T_Optimized']==minT])
						candidates = [sol for sol in hc_sol if data[sol]['T_Optimized']==minT and data[sol]['Npart']==minN]
						cd = data[candidates[0]]
						f_out.write('\t'.join(['5-input', benchmark, cd['Npart'], cd['Sol'], cd['Nodes'], cd['Constraint'], cd['Valid Motif_METIS'], cd['Valid Motif_Optimized'], cd['Cycle Free_METIS'], cd['Cycle Free_Optimized'], cd['T_Metis'], cd['T_Optimized'], cd['cut_Metis'], cd['cut_Optimized']])+'\n') 
					else: 
						# choose from lc sols if no hc sol
						minT = min([data[sol]['T_Optimized'] for sol in lc_sol])
						minN = min([data[sol]['Npart'] for sol in lc_sol if data[sol]['T_Optimized']==minT])
						candidates = [sol for sol in lc_sol if data[sol]['T_Optimized']==minT and data[sol]['Npart']==minN]
						cd = data[candidates[0]]
						f_out.write('\t'.join(['5-input', benchmark, cd['Npart'], cd['Sol'], cd['Nodes'], cd['Constraint'], cd['Valid Motif_METIS'], cd['Valid Motif_Optimized'], cd['Cycle Free_METIS'], cd['Cycle Free_Optimized'], cd['T_Metis'], cd['T_Optimized'], cd['cut_Metis'], cd['cut_Optimized']])+'\n') 
				
				else: 
					print(benchmark)

def load_timestep (filename):
	lines = [open(filename, 'r').read().strip("\n")][0].split('\n')
	timeList = []
	for line in lines[1:]:
		times = line.split('\t')[1].split(',')
		timeList = timeList + times
	return list(filter(None, timeList))

def load_median_connectivity (input_n):
	bm_path = PATH + 'runs/benchmark/'+input_n+'-input-boolean-circuits/'
	rs_path = PATH + 'runs/results/'+input_n+'-input-boolean-circuits/'
	connDict = {}
	for benchmark in os.listdir(rs_path):
		if benchmark.isdigit():	
			edgelist = bm_path + benchmark + '/DAG.edgelist'
			G = load_graph (edgelist)
			in_nodes, out_nodes, nonprimitives  = gp.get_nonprimitive_nodes (G)
			G_primitive = gp.get_G_primitive (G, nonprimitives)
			connDict[benchmark] = {}
			for npart in os.listdir(rs_path + benchmark + '/nparts/'):
				if npart.isdigit():
					part_sol = rs_path + benchmark + '/nparts/' + npart + '/part_solns.txt'
					cut, partDict = gp.load_metis_part_sol (part_sol)
					part_opt = (cut, [gp.get_part(partDict, n) for n in G_primitive.nodes()])
					matrix, partG = gp.partition_matrix (G_primitive, part_opt)
					# calculate the median connectivity
					sum_row = matrix.sum(axis=1)
					sum_col = matrix.sum(axis=0)
					mean_degree = np.median(sum_col + sum_row.T)
					connDict[benchmark][npart] = mean_degree
	return connDict


def plot_timestep ():
	""" plot the distribution of timesteps at which improvement is made """
	bm_path = PATH + 'runs/benchmark/5-input-boolean-circuits/'
	rs_path = PATH + 'runs/results/5-input-boolean-circuits/'
	hc_timeDict, lc_timeDict = {}, {}
	hc_timeList, lc_timeList = [], []
	connectivityDict = load_median_connectivity ('5')
	hc_avgConn_opt_List, lc_avgConn_opt_List = [], []   # associate average degree of connectivity among subgraphs, and number of optimized steps
	hc_avgNode_opt_List, lc_avgNode_opt_List = [], []   # associate number of nodes in partitioned subgraphs, and number of optimized steps
	for benchmark in os.listdir(rs_path):
		if benchmark.isdigit():	
			for npart in os.listdir(rs_path + benchmark + '/nparts/'):
				if npart.isdigit():
					hc_opt_file = rs_path+benchmark+'/nparts/'+npart+'/optimized_hc/'+'part improved.txt'
					lc_opt_file = rs_path+benchmark+'/nparts/'+npart+'/optimized_lc/'+'part improved.txt'
					G = load_graph (bm_path + benchmark + '/DAG.edgelist')
					primN = len(list(G.nodes())) - 6 # number of primitive vertices
					if os.path.exists(hc_opt_file):
						timeList = load_timestep (hc_opt_file)
						hc_timeList = hc_timeList + timeList 
						if npart in hc_timeDict:
							hc_timeDict[npart] = hc_timeDict[npart] + timeList
						else: 
							hc_timeDict[npart] = timeList
						hc_avgNode_opt_List.append((primN/int(npart), len(timeList)))
						hc_avgConn_opt_List.append((connectivityDict[benchmark][npart], len(timeList)))
					if os.path.exists(lc_opt_file):
						timeList = load_timestep (lc_opt_file)
						lc_timeList = lc_timeList + timeList 
						if npart in lc_timeDict:
							lc_timeDict[npart] = lc_timeDict[npart] + timeList
						else: 
							lc_timeDict[npart] = timeList
						lc_avgNode_opt_List.append((primN/int(npart), len(timeList)))
						lc_avgConn_opt_List.append((connectivityDict[benchmark][npart], len(timeList)))

	### plot distribution of optimization at timesteps (grouped by constraints, and N)
	# fig = plt.figure(figsize=(9,9)) 
	# # ax.hist([int(t) for t in hc_timeList], fc='orange', histtype='step', label='High Constraint')
	# # ax.hist([int(t) for t in lc_timeList], fc='blue', histtype='step', label='Low Constraint')
	# idx = 1
	# for npart in ['3', '4', '5', '6', '7', '8', '9', '10', '11']:
	# 	ax  = fig.add_subplot(3,3,idx)
	# 	ax.hist([int(t) for t in lc_timeDict[npart]], histtype='step', label = npart)
	# 	ax.set_ylim([0, 2000])
	# 	ax.set_xlim([0, 10000])
	# 	if idx != 7: 
	# 		ax.set_xticks([])
	# 		ax.set_yticks([])
	# 	if idx == 7:
	# 		ax.set_xlabel('Time Step')
	# 		ax.set_ylabel('Count')
	# 	ax.set_title('N='+npart)
	# 	for spine in ['top', 'right']:
	# 		ax.spines[spine].set_linewidth(0)
	# 	idx += 1
	# plt.savefig(PATH+'runs/results/analysis/Optimization at Timestep (lc grouped by N).pdf', dpi=200)
	# plt.show()

	### plot average node in subgraphs vs. number of optimized steps
	fig = plt.figure( figsize=(5,5) ) 
	ax  = fig.add_subplot (111)
	x = list(list(zip(*hc_avgNode_opt_List))[0])
	y = list(list(zip(*hc_avgNode_opt_List))[1])
	ax.plot (x, y, 'o', c='#1f77b4', fillstyle='none', label='High Constraint')
	x = list(list(zip(*lc_avgNode_opt_List))[0])
	y = list(list(zip(*lc_avgNode_opt_List))[1])
	ax.plot (x, y, 'o', c='#ff7f0f', fillstyle='none', label='Low Constraint')
	ax.set_xlim([2.5, 7.5])
	ax.set_ylim([0, 80])
	ax.set_xlabel('Average vertices per subgraph')
	ax.set_ylabel('Number of optimization')
	plt.legend(loc='upper right', frameon=False)
	plt.savefig(PATH+'runs/results/analysis/Number of optimization vs. Average Nodes.pdf', dpi=200)
	plt.show()

	### plot number of optimization vs. median degree of connectivity 
	fig = plt.figure( figsize=(5,5) ) 
	ax  = fig.add_subplot (111)
	x = list(list(zip(*hc_avgConn_opt_List))[0])
	y = list(list(zip(*hc_avgConn_opt_List))[1])
	ax.plot (x, y, 'o', c='#1f77b4', fillstyle='none', label='High Constraint')
	x = list(list(zip(*lc_avgConn_opt_List))[0])
	y = list(list(zip(*lc_avgConn_opt_List))[1])
	ax.plot (x, y, 'o', c='#ff7f0f', fillstyle='none', label='Low Constraint')
	# ax.set_xlim([2.5, 7.5])
	ax.set_ylim([0, 80])
	ax.set_xlabel('Subgraph median degree of connectivity')
	ax.set_ylabel('Number of optimization')
	plt.legend(loc='upper right', frameon=False)
	plt.savefig(PATH+'runs/results/analysis/Number of optimization vs. Subgraph connectivity.pdf', dpi=200)
	plt.show()	
	

def to_tuple (matrix):
	return tuple(tuple(arr.tolist()) for arr in matrix)

def count_motif ():

	motifDict = {}

	for N in [4, 5]:
		bm_path = PATH + 'runs/benchmark/'+str(N)+'-input-boolean-circuits/'
		rs_path = PATH + 'runs/results/'+str(N)+'-input-boolean-circuits/'
		best_sol = PATH + 'runs/results/'+str(N)+'-input-boolean-circuits/best_solns.txt'
		data = load_sol (best_sol)
		for benchmark in data.keys():
			print(benchmark)
			edgelist = bm_path + benchmark + '/DAG.edgelist'
			G = load_graph (edgelist)
			in_nodes, out_nodes, nonprimitives  = gp.get_nonprimitive_nodes (G)
			G_primitive = gp.get_G_primitive (G, nonprimitives)

			if int(data[benchmark]['Sol']) != 0:
				part_sol = rs_path + benchmark + '/nparts/' + data[benchmark]['Npart'] + '/optimized_' + data[benchmark]['Constraint'] + '/part_solns.txt'
				solDict = gp.load_opt_part_sol (part_sol)
				cut = int(solDict[int(data[benchmark]['Sol'])]['cut'])
				part = solDict[int(data[benchmark]['Sol'])]['part']
				part_opt = (cut, [gp.get_part(part, n) for n in G_primitive.nodes()])
				matrix, partG = gp.partition_matrix (G_primitive, part_opt)
			else: 
				part_sol = rs_path + benchmark + '/nparts/' + data[benchmark]['Npart'] + '/part_solns.txt'
				cut, partDict = gp.load_metis_part_sol (part_sol)
				part_opt = (cut, [gp.get_part(partDict, n) for n in G_primitive.nodes()])
				matrix, partG = gp.partition_matrix (G_primitive, part_opt)

			# count the occurrence of each subnetwork
			for cell in list(partG.nodes()):
				neighbors = []
				for e in list(partG.edges()):
					c1, c2 = e[0], e[1]
					if e[0] == cell: neighbors.append(e[1])
					elif e[1] == cell: neighbors.append(e[0])
					subG_cells = list(set(neighbors)) + [cell]
					subG_matrix, subG_partG = gp.get_part_matrix_subG (matrix, partG, subG_cells)

				dim = subG_matrix.shape[0]
				# print('part matrix', matrix)
				I = np.identity(n=dim)  # identity matrix
				dmList = []
				for pm in itertools.permutations(I):  # permutation matrix
					pm = np.array(pm)
					# print('perm matrix', pm)
					dm = np.dot(pm, subG_matrix)               # permutate by rows
					# print('product', np.dot(dm, pm.T))
					dm = np.dot(dm, pm.T)
					dmList.append (to_tuple(dm))              # append all diagnizable matrix of partition matrix 
				# print('all diag matrix', dmList)

				# motif_exist = set(dmList).intersect(set(motifDict.keys()))
				if any(dm in motifDict for dm in dmList) == False:
					# print('motif not found in dict', matrix)
					motifDict[to_tuple(subG_matrix)] = (1/dim)/len(list(data.keys()))
					# print('new motifdict', motifDict)
				else:
					motif = set(dmList).intersection(set(motifDict.keys()))
					# print('motif', motif,'already in dict')
					motifDict[tuple(motif)[0]] += (1/dim)/len(list(data.keys()))


	# write output 
	f_out = open(PATH + 'runs/results/analysis/motif freq.txt', 'w')
	f_out.write('Motif matrix\tOccur\n')
	for motif in motifDict:
		f_out.write(str(motif)+'\t'+str(motifDict[motif])+'\n')
		print(motif, motifDict[motif])


def load_motif_data (inputfile):
 	lines = [open(inputfile, 'r').read().strip("\n")][0].split('\n')
 	motif_dict = {}
 	for line in lines[1:]:
 		tokens = line.split('\t')
 		motif_dict[tokens[0]] = float(tokens[1])
 	return motif_dict

def generate_comm_graph (matrix):
	# generate DAG representing cell-cell communication from the adjency matrix
	rows, cols = np.where(matrix != 0)
	edges = zip(rows.tolist(), cols.tolist())
	partG = nx.DiGraph() 
	partG.add_edges_from(edges) 
	return partG

def plot_motif_occurence (inputfile):
	""" plot the 10 highest occuring motifs from motif frequency results """
	data = load_motif_data (inputfile)
	sort_occur = sorted(data, key=data.get, reverse=True)  # sort occur from highest to lowest
	sort_occur_freq = [data[k]/sum(data.values()) for k in sort_occur]

	## plot motif
	# fig = plt.figure(figsize=(12,12))
	# sp = 1
	# for motif in sort_occur:
	# 	print(sp, motif)
	# 	if sp <= 25:
	# 		ax = fig.add_subplot(5,5,sp)
	# 		matrix = np.array(eval(motif))
	# 		partG = generate_comm_graph(matrix)
	# 		pos = graphviz_layout(partG, prog='dot')
	# 		nx.draw_networkx_nodes(partG, pos, node_size=30, node_color='#b18ea6')
	# 		labels={n:n for n in partG.nodes()}
	# 		nx.draw_networkx_edges(partG, pos)
	# 		ax.axis('off')
	# 	sp += 1

	# plt.savefig(PATH + 'runs/results/analysis/motif occur.pdf', dpi=200)
	# plt.show()

	## plot frequency 
	# fig = plt.figure(figsize=(8,3))
	# ax  = fig.add_subplot(111)
	# x = list(range(1, len(sort_occur_freq)+1))
	# ax.bar(x, sort_occur_freq, color='white', edgecolor='k')
	# plt.axvline(x=25.5, color='r', linestyle='--')
	# ax.set_xlabel('Subgraph Network Motif')
	# ax.set_ylabel('Frequency (%)')
	# plt.legend(loc='upper right', frameon=False)
	# for spine in ['top', 'right']:
	# 	ax.spines[spine].set_linewidth(0)
	# plt.subplots_adjust(bottom=0.15)
	# plt.savefig(PATH + 'runs/results/analysis/motif frequncy.pdf', dpi=200)
	# plt.show()
	# print(sum(sort_occur_freq[0:25]),len(sort_occur_freq[0:25]), sum(sort_occur_freq))


def load_deltaD (N):
	best_soln = PATH + 'runs/results/'+str(N)+'-input-boolean-circuits/best_solns.txt'
	data = load_sol	(best_soln)
	deltaD = []
	for bm in data: 
		D_ori = data[bm]['T_Metis']
		D_opt = data[bm]['T_Optimized']
		deltaD.append (int(D_ori) - int(D_opt))

	return deltaD

def plot_deltaD ():
	""" plot delta T of all 4- and 5-input circuits """

	deltaD_4 = load_deltaD (4)
	deltaD_5 = load_deltaD (5)

	fig = plt.figure() 
	ax  = fig.add_subplot(111)
	colors = ['orange', 'blue']
	labels = ['4-input Boolean circuits', '5-input Boolean circuits']
	plt.hist([deltaD_4, deltaD_5], 10, histtype='step', stacked=False, fill=False, label=labels)
	# ax.hist(count1, ec='k', histtype='step', fill=False, facecolor='wheat')
	# ax.hist(count2, ec='k', histtype='step', fill=False, facecolor='g')
	ax.set_xlabel('Delta Depth')
	ax.set_ylabel('Count')
	plt.legend(loc='upper right', frameon=False)
	plt.gca().yaxis.set_major_locator(mticker.MultipleLocator(4))
	plt.gca().xaxis.set_major_locator(mticker.MultipleLocator(5))
	for spine in ['top', 'right']:
		ax.spines[spine].set_linewidth(0)
	# plt.savefig('Node distribution.pdf', dpi=200)
	plt.show()


def visualize_subnetworks_unmet_constraint (path, constraint):
	""" for each optimization attempt, visualize the nsubnetworks that unmet constraints """
	bm_path = path + 'runs/benchmark/electronic-circuits/alu/'
	sol_path = path + 'runs/results/electronic-circuits/alu/nparts/46/'
	npart = 46
	# for npart in os.listdir(sol_path):
	# 	if npart.startswith('.') == False and npart != 'best_solns.txt':
	# 		print ('npart', npart)
	# load graph and metis partition
	edgelist = bm_path + '/DAG.edgelist'
	G = load_graph (edgelist)
	in_nodes, out_nodes, nonprimitives  = gp.get_nonprimitive_nodes (G)
	G_primitive = gp.get_G_primitive (G, nonprimitives)

	# visualize original partition
	cut, partDict = gp.load_metis_part_sol (sol_path+'part_solns.txt')
	part_opt = [gp.get_part(partDict, n) for n in G_primitive.nodes()]
	matrix, partG = gp.partition_matrix (G_primitive, part_opt)
	cell_unmet_const, cell_met_const = gp.get_cells_unmet_constraint (matrix, partG, [5], 'TRUE')
	gp.visualize_assignment_graphviz (G, part_opt, nonprimitives, 'TRUE', sol_path, 0, cell_unmet_const)


	# part_sol = sol_path + npart + '/part_solns.txt'
	part_sol = sol_path + 'part_solns.txt'
	cut, partDict = gp.load_metis_part_sol (part_sol)

	# opt_file = sol_path + npart + '/optimized_'+constraint+'/part_solns.txt'
	# opt_sol_file = sol_path + npart + '/optimized_'+constraint+'/part_solns.txt'
	opt_file = sol_path + 'optimized_lc/part_solns.txt'
	if os.path.exists (opt_file):
		print('opt exists')
		timeList = load_timestep (opt_file)
		if timeList != []:
			solDict = gp.load_opt_part_sol (opt_file)
			for iteration in solDict.keys():
				part = solDict[iteration]['part']
				if part != partDict:
					part_opt = [gp.get_part(part, n) for n in G_primitive.nodes()]
					matrix, partG = gp.partition_matrix (G_primitive, part_opt)
					cell_unmet_const, cell_met_const = gp.get_cells_unmet_constraint (matrix, partG, [5], 'TRUE')
					print(iteration, solDict[iteration]['T'], len(cell_unmet_const)/int(npart))
					# gp.visualize_assignment_graphviz (G, part_opt, nonprimitives, 'TRUE', sol_path+'/optimized_'+constraint, iteration, cell_unmet_const)


def generate_edgelist_of_cells (outdir): 
	""" for partitioned cells, generate an edgelist """
	bm_path = outdir + 'runs/benchmark/electronic-circuits/md5Core/'
	sol_path = outdir + 'runs/results/electronic-circuits/md5Core/nparts/40-minimizing unmet/optimized_lc/'
	part_sol = sol_path + 'part_solns.txt'

	edgelist = bm_path + '/DAG.edgelist'
	G = load_graph (edgelist)
	in_nodes, out_nodes, nonprimitives  = gp.get_nonprimitive_nodes (G)
	G_primitive = gp.get_G_primitive (G, nonprimitives)

	solDict = gp.load_opt_part_sol (part_sol)
	iteration = 11

	part = solDict[iteration]['part']
	part_opt = [gp.get_part(part, n) for n in G_primitive.nodes()]

	node_to_partDict = {}
	for n in G_primitive.nodes():
		node_to_partDict[n] = str(gp.get_part(part, n))

	matrix, partG = gp.partition_matrix (G_primitive, part_opt)

	# nx.write_edgelist (partG, sol_path+'md5Core_part40_iteration11.edgelist')

	f_out = open (sol_path + 'neighbors.txt', 'w')
	for node in G_primitive.nodes():
		neighbors = G_primitive.neighbors(node)
		f_out.write(node + '\t' + ','.join(neighbors)+'\n')

	return node_to_partDict

def get_cut_edges (g, node_to_partDict):

	cut_edges, internal_edges = [], []
	for e in g.edges():
		n1, n2 = e[0], e[1]
		if node_to_partDict[n1] != node_to_partDict[n2]: 
			cut_edges.append(e)
		else: 
			internal_edges.append(e)
	return cut_edges, internal_edges

def _scale_xy (array, x_scaler, y_scaler):
	new_array = np.array ([array[0] * x_scaler, array[1] * y_scaler])
	return new_array 


def plot_cell_edgelist(
        input_edgelist_fp: str,
        part_edgelist_fp: str,
        partition: dict, 
        save_file: bool = False,
        output_filename: str = 'test.jpg',
):
	original_network = False
	g = nx.read_edgelist(input_edgelist_fp, create_using=nx.DiGraph())
	in_nodes, out_nodes, nonprimitives  = gp.get_nonprimitive_nodes (g)
	g = gp.get_G_primitive (g, nonprimitives)
	pg = nx.read_edgelist(part_edgelist_fp)
	cut_edges, internal_edges = get_cut_edges (g, partition)

	plt.figure(num=None, figsize=(15, 15), dpi=80)
	img_path = "/Users/jgzhang/Work/Densmore_lab/Partition/code_version/v2/genetic-circuit-partitioning/2021.4/runs/results/Logic-gate-nor-us.png"
	graph_path = "/Users/jgzhang/Work/Densmore_lab/Partition/code_version/v2/genetic-circuit-partitioning/2021.4/runs/results/electronic-circuits//"
	img_list = []
	nor_img = mpimg.imread(img_path)
	for _ in pg.nodes():
		img_list.append(nor_img)

	positions = nx.nx_agraph.graphviz_layout(pg, prog='dot')
	positions = nx.drawing.layout.rescale_layout_dict (positions, scale=10)
	positions = {key: _scale_xy(value, 1,1) for (key, value) in positions.items()}

    ## position of nodes
	pos_communities = dict()
	for node, part in partition.items():
		pos_communities[node] = positions[part]

	# position nodes within partition
	partDict = dict()
	for node, part in partition.items():
		try:
			partDict[part] += [node]
		except KeyError:
			partDict[part] = [node]

	pos_nodes = dict()
	for ci, nodes in partDict.items():
		subgraph = g.subgraph(nodes)
		pos_subgraph = nx.planar_layout(subgraph)
		# pos_subgraph = nx.nx_agraph.graphviz_layout(subgraph)
		pos_nodes.update(pos_subgraph)

	# combine position
	pos = dict()
	for node in g.nodes():
		print('community position', pos_communities[node])
		print('node position', pos_nodes[node])
		print('new position', np.array(pos_communities[node]) + np.array(pos_nodes[node]))
		pos[node] = np.array(pos_communities[node]) + _scale_xy (pos_nodes[node], 0.5, 0.5)

    # Calculate the number of ranks in here so you can figure out how many
    # colors you need...
	y_pos = sorted(list({position[1] for position in positions.values()}))
	sns.color_palette('Set2', len(y_pos))
	colors = [y_pos.index(position[1]) for position in positions.values()]
	plt.title('RCA4 Boolean Logic Network', fontsize=30, ha='center')
	if original_network:
		nx.draw (
			pg, 
			pos=positions,
			with_labels=True,
			node_color=colors,
			width=0,
			node_size=1000,
			node_shape='s',
			linewidths=30,
			horizontalalignment='left',
			verticalalignment='bottom',
			alpha=0.2,
		)
		nx.draw_networkx_edges (g, pos=pos, edgelist=cut_edges, edge_color='red')
		nx.draw_networkx_edges (g, pos=pos, edgelist=internal_edges, edge_color='k')
		nx.draw(
			g,
			pos=pos,
			with_labels=True,
			width=0,
			horizontalalignment='left',
			verticalalignment='bottom',
			alpha=0.7,
		)

		plt.draw()
		plt.savefig(graph_path + 'original_network.jpg')
		return
    # I'm going to rotate the graph so it makes more sense in terms of an
    # electrical circuit.
	flip_xy_pos = {}
	flipped_pos = {node: (-y, x) for (node, (x, y)) in pos.items()}
	flipped_positions = {node: (-y, x) for (node, (x, y)) in positions.items()}
	nx.draw (
		pg, 
		pos=flipped_positions,
		with_labels=False,
		node_color=colors,
		width=0,
		node_size=2000,
		node_shape='s',
		linewidths=40,
		horizontalalignment='left',
		verticalalignment='bottom',
		alpha=0.2,
		)
	nx.draw_networkx_edges (g, pos=flipped_pos, edgelist=cut_edges, edge_color='red')
	nx.draw_networkx_edges (g, pos=flipped_pos, edgelist=internal_edges, edge_color='k')
	nx.draw(
		g,
		pos=flipped_pos,
		with_labels=True,
		width=0,
		horizontalalignment='left',
		verticalalignment='bottom',
		alpha=0.7,
	)


    # position nodes within each partition 


    # plt.imshow(nor_img)
	# ax = plt.gca()
	# fig = plt.gcf()
	# trans = ax.transData.transform
	# trans2 = fig.transFigure.inverted().transform
	# imsize = 0.05  # this is the image size
	# for index, n in enumerate(pg.nodes()):
	# 	(x, y) = flipped_pos[n]
	# 	xx, yy = trans((x, y))  # figure coordinates
	# 	xa, ya = trans2((xx, yy))  # axes coordinates
	# 	print(f'{xa=}')
	# 	print(f'{ya=}')
	# 	print(f'{imsize=}')
	# 	print(f'{xa - imsize / 2.0=}')
	# 	print(f'{ya - imsize / 2.0=}')
	# 	a = plt.axes([(xa / 0.7975) - 0.17, (ya / 0.805) - 0.155, imsize, imsize])
	# 	a.imshow(img_list[index])
	# 	a.set_aspect('equal')
	# 	a.axis('off')
	# plt.show()
	plt.draw()
	plt.savefig(graph_path + 'flipped_network_new.pdf', dpi=300)



if __name__ == '__main__':
	PATH = '/Users/jgzhang/Work/Densmore_lab/Partition/code_version/v2/genetic-circuit-partitioning/2021.4/'
	# count_nodes (PATH + 'runs/results/4-input-boolean-circuits', PATH + 'runs/results/5-input-boolean-circuits')
	# partition_stats ()
	# compile_best_solutions ()
	# plot_timestep ()
	# count_motif()
	# plot_motif_occurence (PATH + 'runs/results/analysis/motif freq.txt')
	# plot_deltaD ()
	# G = load_graph (PATH + 'runs/benchmark/electronic-circuits/hexDisplay/DAG.edgelist')
	# plot_network (G, PATH + 'runs/benchmark/electronic-circuits/hexDisplay/DAG.pdf')

	# best_soln = PATH + 'runs/results/5-input-boolean-circuits/best_solns.txt'
	# data = load_sol	(best_soln)
	# avgNode = np.mean([(int(data[bm]['Nodes'])-6)/int(data[bm]['Npart']) for bm in data])
	# print(avgNode)
	partition = generate_edgelist_of_cells (PATH)
	# visualize_subnetworks_unmet_constraint (PATH, 'lc')

	# plot_cell_edgelist (PATH+'runs/benchmark/electronic-circuits/md5Core/DAG.edgelist', PATH+'runs/results/electronic-circuits/md5Core/nparts/40-minimizing unmet/optimized_lc/md5Core_part40_iteration11.edgelist', partition)

