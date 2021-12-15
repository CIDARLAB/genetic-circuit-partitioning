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
from collections import Counter
from scipy import stats

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

	# pos = nx.kamada_kawai_layout(G)
	# # pos = nx.random_layout(G)
	# plt.figure(num=None, figsize=(3,3), dpi=80)
	# nx.draw(
	#     G,
	#     pos=pos,
	#     horizontalalignment='left',
	#     verticalalignment='bottom',
	#     node_color='powderblue'
	# )
	# # plt.savefig(outfig, dpi=200)
	# plt.show()

	### graphviz layout
	pos = nx.kamada_kawai_layout(G)
	print(len(G.nodes()))
	# pos =  nx.random_layout(G)
	# pos = graphviz_layout(G, prog='dot')
	plt.figure(num=None, figsize=(3,3), dpi=80)
	nx.draw (G, pos, node_color='powderblue', width=0.25, node_size=30, arrowsize=6, with_labels=False)
	plt.savefig (outfig, dpi=200)
	plt.show()

	### plot flipped positions
	# print(len(G.nodes()))
	# plt.figure(num=None, figsize=(4,3), dpi=80)
	# flip_xy_pos = {}
	# flipped_pos = {node: (-y, x) for (node, (x, y)) in pos.items()}
	# nx.draw (G, pos=flipped_pos, node_color='powderblue', width=0.25, node_size=30, arrowsize=6,horizontalalignment='left',verticalalignment='bottom')
	# plt.savefig (outfig, dpi=200)
	# plt.show()

	### plot circular layout
	# plt.figure(figsize=(4,4))
	# pos = nx.circular_layout(G) 
	# nx.draw(G, pos=pos, node_size=30, node_color='powderblue', width=0.25, arrowsize=6, connectionstyle='arc3, rad=0.1')
	# # nx.draw_networkx_nodes(G, pos, node_size=20, node_color='powderblue')
	# # nx.draw_networkx_edges(G, pos, width=0.25, arrowsize=6, connectionstyle='arc3, rad=0.1')
	# plt.savefig(outfig, dpi=200)
	# plt.show()

def plot_centrality (G, outfig):
	""" compute and plot closeness centrality of nodes, betweenness centrality of nodes, and degree of nodes """
	closenessDict = nx.algorithms.centrality.closeness_centrality (G)
	sorted_closeness = sorted(closenessDict.items(), key=lambda x: x[1])[::-1]

	betweennessDict = nx.algorithms.centrality.betweenness_centrality (G)
	sorted_betweenness = sorted(betweennessDict.items(), key=lambda x: x[1])[::-1]
	
	degreeDict = nx.algorithms.centrality.degree_centrality (G) 
	sorted_degree = sorted(degreeDict.items(), key=lambda x: x[1])[::-1]

	fig = plt.figure(figsize=(5,3))
	nodes = list(G.nodes())
	y1 = [closenessDict[x] for x in nodes]
	y1_norm = [(y-min(y1))/(max(y1)-min(y1)) for y in y1]
	zscore1 = stats.zscore(y1)
	y2 = [betweennessDict[x] for x in nodes]
	y2_norm = [(y-min(y2))/(max(y2)-min(y2)) for y in y2]
	zscore2 = stats.zscore(y2)
	y3 = [degreeDict[x] for x in nodes]
	y3_norm = [(y-min(y3))/(max(y3)-min(y3)) for y in y3]
	zscore3 = stats.zscore(y3)
	x = range(1, len(nodes)+1)
	ax1 = fig.add_subplot(311)
	ax1.plot(x, zscore1, 'o-', markersize=4, c='k')
	ax2 = fig.add_subplot(312)
	ax2.plot(x, zscore2, 'o-', markersize=4, c='k')
	ax3 = fig.add_subplot(313)
	ax3.plot(x, zscore3, 'o-', markersize=4, c='k')
	# ax1.set_ylim([0, 0.4])
	ax1.set_xticks([])
	# ax2.set_ylim([0, 0.4])
	ax2.set_xticks([])
	# ax3.set_ylim([0, 0.4])
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
	bm_path = path + 'runs/benchmark/bionetwork/random_graph/n30_p0.04/'
	sol_path = path + 'runs/results/bionetwork/RG_n30_p0.04/nparts/'
	npart = 5
	# for npart in os.listdir(sol_path):
	# 	if npart.startswith('.') == False and npart != 'best_solns.txt':
	# 		print ('npart', npart)
	# load graph and metis partition
	edgelist = bm_path + 'DAG.edgelist'
	G = load_graph (edgelist)
	in_nodes, out_nodes, nonprimitives  = gp.get_nonprimitive_nodes (G)
	# G_primitive = gp.get_G_primitive (G, nonprimitives)
	G_primitive = G

	# visualize original partition
	cut, partDict = gp.load_metis_part_sol (sol_path+str(npart)+'/part_solns.txt')
	part_opt = [gp.get_part(partDict, n) for n in G_primitive.nodes()]
	matrix, partG = gp.partition_matrix (G_primitive, part_opt)
	# qs_matrix = gp.calc_qs (G_primitive, part_opt)
	cell_unmet_const, cell_met_const = gp.get_cells_unmet_constraint (matrix, partG, [1], 'FALSE')
	# cell_unmet_const, cell_met_const = gp.get_cells_unmet_qs_constraint (matrix, partG, qs_matrix, [4], 'TRUE')
	gp.visualize_assignment_graphviz (G, part_opt, nonprimitives, 'FALSE', sol_path, 0, [])

	for conn in os.listdir(sol_path+str(npart)+'/optimized'):
		print(conn)
		if conn == '4':
			part_sol = sol_path+str(npart)+'/part_solns.txt'
			# part_sol = sol_path + 'part_solns.txt'
			cut, partDict = gp.load_metis_part_sol (part_sol)

			# f_out = open (sol_path + 'optimized_lc/iteration-2.txt', 'w')
			opt_file = sol_path+str(npart)+'/optimized/'+conn+'/part_solns.txt'

			solDict = gp.load_opt_part_sol (opt_file)
			for iteration in solDict.keys():
				print('iteration', iteration)
				part = solDict[iteration]['part']
				if part != partDict:
					part_opt = [gp.get_part(part, n) for n in G_primitive.nodes()]
					matrix, partG = gp.partition_matrix (G_primitive, part_opt)
					# qs_matrix = gp.calc_qs (G_primitive, part_opt)
					# median_qs_best = np.mean(np.sum(qs_matrix, axis=1))
					# cell_unmet_const, cell_met_const = gp.get_cells_unmet_qs_constraint (matrix, partG, qs_matrix, [4], 'TRUE')
					cell_unmet_const, cell_met_const = gp.get_cells_unmet_constraint (matrix, partG, [int(conn)], 'FALSE')
					if len(cell_unmet_const) == 0: print ('solution')
					# for idx, p in enumerate(part):
					# 	# print('Partition '+str(part_num)+' '+','.join(part[p]))
					# 	qs = list(qs_matrix[idx])
					# 	sumqs = sum(qs)
					# 	f_out.write('Partition '+str(idx+1)+'\t'+str(sumqs)+'\t'+str(len(part[p]))+'\t'+', '.join([str(int(v)) for v in qs])+'\t'+', '.join(part[p])+'\t'+'\t'+', '.join([v for v in part[p]])+'\n')
					# 	print('Partition '+str(idx+1)+'\t'+str(sumqs)+'\t'+str(len(part[p]))+'\t'+', '.join([str(int(v)) for v in qs])+'\t'+', '.join(part[p])+'\t'+'\t'+', '.join([v for v in part[p]])+'\n')
					# print(iteration, median_qs_best, solDict[iteration]['T'], len(cell_unmet_const))
					gp.visualize_assignment_graphviz (G, part_opt, nonprimitives, 'FALSE', sol_path+str(npart)+'/optimized/'+conn, iteration, cell_unmet_const)



def compare_runs (path, constraint):
	""" for each optimization attempt, visualize the nsubnetworks that unmet constraints """
	bm_path = path + 'runs/benchmark/electronic-circuits/alu/'

	# for npart in os.listdir(sol_path):
	# 	if npart.startswith('.') == False and npart != 'best_solns.txt':
	# 		print ('npart', npart)
	# load graph and metis partition
	edgelist = bm_path + '/DAG.edgelist'
	G = load_graph (edgelist)
	in_nodes, out_nodes, nonprimitives  = gp.get_nonprimitive_nodes (G)
	G_primitive = gp.get_G_primitive (G, nonprimitives)

	nparts = [47, 48, 49, 50]

	dataDict = {}

	f_out = open (path + 'runs/results/electronic-circuits/alu/nparts/Compare Runs Results.txt', 'w')
	f_out.write('Cell number\tIteration\tEnd Cell Number\tMedian QS Per Cell\tNumber of Cells with >=3 QS\tFraction of Cells with >=3 QS\tNetwork Depth\n')
	# plot constrant unmet vs. 
	for idx, npart in enumerate(nparts):
		print(npart)
		sol_path = path + 'runs/results/electronic-circuits/md5Core/nparts/'+str(npart)
		part_sol = sol_path + '/part_solns.txt'
		cut, partDict = gp.load_metis_part_sol (part_sol)
		opt_file = sol_path + '/optimized_lc/part_solns.txt'
		if os.path.exists (opt_file):
			x, y = [], []
			timeList = load_timestep (opt_file)
			if timeList != []:
				solDict = gp.load_opt_part_sol (opt_file)
				for iteration in solDict.keys():
					part = solDict[iteration]['part']
					if part != partDict:
						part_opt = [gp.get_part(part, n) for n in G_primitive.nodes()]
						endN = len(part.keys())
						matrix, partG = gp.partition_matrix (G_primitive, part_opt)
						qs_matrix = gp.calc_qs (G_primitive, part_opt)
						median_qs_best = np.mean([x for x in np.sum(qs_matrix, axis=1) if x>0])
						cell_unmet_const, cell_met_const = gp.get_cells_unmet_qs_constraint (matrix, partG, qs_matrix, [3], 'TRUE')
						# print(iteration, median_qs_best, solDict[iteration]['T'], len(cell_unmet_const), len(cell_unmet_const)/endN)
						f_out.write('\t'.join([str(npart), str(iteration), str(endN), str(round(median_qs_best, 2)), str(len(cell_unmet_const)), str(round(len(cell_unmet_const)/endN, 2)), str(solDict[iteration]['T'])])+'\n')
					
						if endN not in dataDict: 
							dataDict[endN] = {}
							dataDict[endN][1] = {'qs': median_qs_best, 'unmet': len(cell_unmet_const), 'unmetf': len(cell_unmet_const)/endN}
						else: 
							dataDict[endN][max(dataDict[endN].keys())+1] = {'qs': median_qs_best, 'unmet': len(cell_unmet_const), 'unmetf': len(cell_unmet_const)/endN}

	colors = sns.color_palette("husl", len(dataDict.keys()))

	fig = plt.figure (figsize=(7,5))
	ax = fig.add_subplot(111)

	for idx, N in enumerate(sorted(dataDict)): 
		print(idx, N)
		x, y = [], []
		c = colors[idx]
		for itr in sorted(dataDict[N]): 
			x.append (dataDict[N][itr]['qs'])
			y.append (dataDict[N][itr]['unmetf'])
		ax.plot (x, y, marker='o', markersize=3, linestyle='', c=c, label=N)

	ax.set_xlabel('Average QS per cell')
	ax.set_ylabel('Fraction of Cells with >= 3 QS')

	plt.legend(bbox_to_anchor=(1.05, 1))
	fig.subplots_adjust (right=0.7)
	plt.savefig(path + 'runs/results/electronic-circuits/md5Core/nparts/Compare Runs QS vs Fraction of cells with >=3QS.pdf', dpi=200)
	plt.show()





def get_gatenum_distribution (partDict):
	gatenumDict = {}
	gates = list(partDict.values())
	gatenum = [len(v) for v in gates]
	min_gatenum, max_gatenum = min(gatenum), max(gatenum)
	for n in range(min_gatenum, max_gatenum+1):
		count = gatenum.count(n)
		gatenumDict[n] = count
	return gatenumDict


def compare_gatenum_distribution (path, constraint):
	""" 
	compare the gate number distribution in each run 
	"""
	bm_path = path + 'runs/benchmark/electronic-circuits/md5Core/'

	# for npart in os.listdir(sol_path):
	# 	if npart.startswith('.') == False and npart != 'best_solns.txt':
	# 		print ('npart', npart)
	# load graph and metis partition
	edgelist = bm_path + '/DAG.edgelist'
	G = load_graph (edgelist)
	in_nodes, out_nodes, nonprimitives  = gp.get_nonprimitive_nodes (G)
	G_primitive = gp.get_G_primitive (G, nonprimitives)

	nparts = [47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60]

	GateNumDict = {}

	f_out = open (path + 'runs/results/electronic-circuits/md5Core/nparts-run2/Compare Runs Gate Number Distribution.txt', 'w')
	f_out.write('Start cell number\tIteration\tEnd cell number\tGate Number Per Cell\tGate Number Distribution\n')
	# plot constrant unmet vs. 
	for idx, npart in enumerate(nparts):
		print(npart)
		sol_path = path + 'runs/results/electronic-circuits/md5Core/nparts-run2/'+str(npart)
		part_sol = sol_path + '/part_solns.txt'
		cut, partDict = gp.load_metis_part_sol (part_sol)
		opt_file = sol_path + '/optimized_lc/part_solns.txt'
		if os.path.exists (opt_file):
			timeList = load_timestep (opt_file)
			if timeList != []:
				solDict = gp.load_opt_part_sol (opt_file)
				for i_idx, iteration in enumerate(solDict.keys()):
					part = solDict[iteration]['part']
					endN = len(part.keys())
					gatecounts = get_gatenum_distribution (part)
					x = list(gatecounts.keys())
					y = list(gatecounts.values())
					f_out.write('\t'.join([str(npart), str(iteration), str(endN), str([str(v for v in x)]), str([str(v for v in y)])])+'\n')
					if endN not in GateNumDict: 
						GateNumDict[endN] = {}
						GateNumDict[endN][1] = gatecounts
					else: 
						GateNumDict[endN][max(GateNumDict[endN].keys())+1] = gatecounts 

	maxIter = max([max(GateNumDict[N].keys()) for N in GateNumDict])
	colors = sns.color_palette("husl", maxIter)

	fig = plt.figure (figsize=(10,3))

	for idx, N in enumerate(sorted(GateNumDict), 1): 
		print(idx, N)
		ax = fig.add_subplot(2,6,idx)
		for i, itr in enumerate(sorted(GateNumDict[N])): 
			x = GateNumDict[N][itr].keys()
			y = normalize_data( GateNumDict[N][itr].values() )
			c = colors[i]
			ax.plot (x, y, marker='o', markersize=3, c=c, label=itr)
		ax.set_title(str(N))
		ax.set_xlim([0, 6])
		ax.set_ylim([-0.1, 1.1])
		if idx != 7: 
			ax.set_xticks([])
			ax.set_yticks([])
		else: 
			ax.set_xticks([1,2,3,4,5])

	fig.subplots_adjust (hspace=0.3)
	plt.savefig(path + 'runs/results/electronic-circuits/md5Core/nparts-run2/Compare Runs Gate Number Distribution.pdf', dpi=200)
	plt.show()



def get_qs_distribution (qs_matrix):
	qsDist = {}
	qs_sum = list(np.sum(qs_matrix, axis=1))
	for e in qs_sum: 
		if e != 0:
			if e not in qsDist: 
				qsDist[e] = 1
			else: 
				qsDist[e] += 1
	return dict(sorted(qsDist.items()))

def normalize_data (data):
	""" normalize data to within 0-1 region """
	return [(d-min(data))/(max(data)- min(data))for d in data]

def compare_qs_distribution (path, constraint):
	""" 
	compare the qs distribution in each run 
	"""
	bm_path = path + 'runs/benchmark/electronic-circuits/alu/'

	# for npart in os.listdir(sol_path):
	# 	if npart.startswith('.') == False and npart != 'best_solns.txt':
	# 		print ('npart', npart)
	# load graph and metis partition
	edgelist = bm_path + '/DAG.edgelist'
	G = load_graph (edgelist)
	in_nodes, out_nodes, nonprimitives  = gp.get_nonprimitive_nodes (G)
	G_primitive = gp.get_G_primitive (G, nonprimitives)

	nparts = [45]

	QSDict = {}

	f_out = open (path + 'runs/results/electronic-circuits/alu/nparts/Compare QS Distribution.txt', 'w')
	f_out.write('Start cell number\tIteration\tEnd cell number\tQS Number\tQS Distribution by cell\n')
	# plot constrant unmet vs. 
	for idx, npart in enumerate(nparts):
		print(npart)
		sol_path = path + 'runs/results/electronic-circuits/alu/nparts/'+str(npart)
		part_sol = sol_path + '/part_solns.txt'
		cut, partDict = gp.load_metis_part_sol (part_sol)
		opt_file = sol_path + '/optimized_lc/part_solns.txt'
		if os.path.exists (opt_file):
			timeList = load_timestep (opt_file)
			if timeList != []:
				solDict = gp.load_opt_part_sol (opt_file)
				for i_idx, iteration in enumerate(solDict.keys()):
					part = solDict[iteration]['part']
					endN = len(part.keys())
					part_opt = [gp.get_part(part, n) for n in G_primitive.nodes()]
					matrix, partG = gp.partition_matrix (G_primitive, part_opt)
					qs_matrix = gp.calc_qs (G_primitive, part_opt)
					qs_dist = get_qs_distribution (qs_matrix)
					print(list(qs_dist.keys()), list(qs_dist.values()))
					f_out.write('\t'.join([str(npart), str(iteration), str(len(part.keys())), ','.join([str(v) for v in list(qs_dist.keys())]), ','.join([str(v) for v in list(qs_dist.values())])])+'\n')
					if endN not in QSDict: 
						QSDict[endN] = {}
						QSDict[endN][1] = qs_dist 
					else: 
						QSDict[endN][max(QSDict[endN].keys())+1] = qs_dist 
	
	# for each ##endN##, plot QS distribution 
	maxIter = max([max(QSDict[N].keys()) for N in QSDict])
	colors = sns.color_palette("husl", maxIter)

	fig = plt.figure (figsize=(10,3))

	for idx, N in enumerate(sorted(QSDict), 1): 
		print(idx, N)
		ax = fig.add_subplot(2,6,idx)
		for i, itr in enumerate(sorted(QSDict[N])): 
			x = QSDict[N][itr].keys()
			y = normalize_data( QSDict[N][itr].values() )
			c = colors[i]
			ax.plot (x, y, marker='o', markersize=3, c=c, label=itr)
		ax.set_title(str(N))
		ax.set_xlim([0, 8])
		ax.set_ylim([-0.1, 1.1])
		if idx != 7: 
			ax.set_xticks([])
			ax.set_yticks([])
		else: 
			ax.set_xticks([1,2,3,4,5,6,7,8])

	fig.subplots_adjust (hspace=0.3)
	plt.savefig(path + 'runs/results/electronic-circuits/alu/nparts/Compare Runs QS Distribution.pdf', dpi=200)
	plt.show()


def generate_edgelist_of_cells (outdir): 
	""" for partitioned cells, generate an edgelist """
	bm_path = outdir + 'runs/benchmark/electronic-circuits/alu/'
	sol_path = outdir + 'runs/results/electronic-circuits/alu/nparts/45/optimized_lc/'
	part_sol = sol_path + 'part_solns.txt'

	edgelist = bm_path + '/DAG.edgelist'
	G = load_graph (edgelist)
	in_nodes, out_nodes, nonprimitives  = gp.get_nonprimitive_nodes (G)
	G_primitive = gp.get_G_primitive (G, nonprimitives)

	solDict = gp.load_opt_part_sol (part_sol)
	iteration = 2

	part = solDict[iteration]['part']
	part_opt = [gp.get_part(part, n) for n in G_primitive.nodes()]

	node_to_partDict = {}
	for n in G_primitive.nodes():
		node_to_partDict[n] = str(gp.get_part(part, n))

	matrix, partG = gp.partition_matrix (G_primitive, part_opt)

	# nx.write_edgelist (partG, sol_path+'part_iteration2.edgelist')

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
	graph_path = "/Users/jgzhang/Work/Densmore_lab/Partition/code_version/v2/genetic-circuit-partitioning/2021.4/runs/results/electronic-circuits/alu/"
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
	# plt.title('RCA4 Boolean Logic Network', fontsize=30, ha='center')
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
		with_labels=True,
		node_color=colors,
		width=0,
		node_size=2000,
		node_shape='s',
		linewidths=40,
		horizontalalignment='left',
		verticalalignment='bottom',
		alpha=0.2,
		)
	nx.draw_networkx_edges (g, pos=flipped_pos, edgelist=cut_edges, edge_color='k')
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
	plt.savefig(graph_path + 'flipped_network_with_label.pdf', dpi=300)


def convert_name (inputfile):
	""" convert MD5 gate name to Jai's names"""
	lines = [open(inputfile, 'r').read().strip("\n")][0].split('\n')
	converted_names = {}
	for line in lines: 
		tokens = line.split('\t')
		converted_names[tokens[0]] = tokens[1]
	return converted_names


def parameter_scan_bionetworks (PATH):
	""" plot solutions scanned by initial n and connectivity constraint """
	bm = PATH + 'runs/benchmark/bionetwork/random_graph/n30_p0.05/DAG.edgelist'
	sol_path = PATH + 'runs/results/bionetwork/RG_n30_p0.05/nparts/'

	G = load_graph (bm)
	in_nodes, out_nodes, nonprimitives  = gp.get_nonprimitive_nodes (G)
	G_primitive = gp.get_G_primitive (G, nonprimitives)
	print(nonprimitives)

	sols = []
	# for targetn in os.listdir(sol_path):
	for targetn in ['8']:
		if targetn.startswith('.') == False and targetn.isdigit():
			print('target n', targetn)
			part_sol = sol_path + targetn + '/part_solns.txt'
			cut, partDict = gp.load_metis_part_sol (part_sol)
			opt_path = sol_path + targetn + '/optimized/'
			if os.path.exists(opt_path):
				# for conn in os.listdir(opt_path):
				for conn in ['7']:
					if conn.startswith('.')  == False and conn.isnumeric(): 
						print('connectivity', conn)
						opt_part_sol = opt_path + conn + '/part_solns.txt'
						if os.path.exists(opt_part_sol):
							solDict = gp.load_opt_part_sol (opt_part_sol)
							for ite in solDict.keys():
								print(solDict[ite])
								part = solDict[ite]['part']
								part_opt = (cut, [gp.get_part(part, n) for n in G.nodes()])
								for n in list(G.nodes()):
									print(n, gp.get_part(part, n))
								endN = len(set(part_opt[1]))
								matrix, partG = gp.partition_matrix (G, part_opt[1])
								motif_allowed = gp.check_motif_allowed(matrix, conn)
								if motif_allowed: 
									print('iteration', ite)
									sols.append((endN, int(conn)))
									gp.visualize_assignment_graphviz (G, part_opt[1], nonprimitives, 'FALSE', opt_path + str(conn) + '/', 'iter_'+str(ite)+'_N_'+str(endN)+'_conn_'+str(conn), [])
	
	# plot solutions
	# fig = plt.figure(figsize=(5,5))
	# ax = fig.add_subplot(111)
	# sol_count = []
	# counted = []
	# for sol in sols: 
	# 	if sol not in counted: 
	# 		count = sols.count(sol)
	# 		sol_count.append ((sol, count))
	# 		counted.append(sol)
	# print(sol_count)
	# x, y, z = [], [], []
	# for sol in sol_count: 
	# 	x.append (sol[0][0])
	# 	y.append (sol[0][1])
	# 	z.append (sol[1])
	# print(x)
	# print(y)
	# print(z)
	# ax.scatter(np.array(x), np.array(y), s=np.array(z)*20, alpha=0.4, c='blue', edgecolors='grey', linewidth=1)
	# # ax.set_xlim([0, 10])
	# # ax.set_ylim([0, 7])
	# ax.set_xlabel('Number of Submodules in Valid Solutions')
	# ax.set_ylabel('Maximum Connuections between submodules')
	# plt.savefig(sol_path+'Solution_space.png', dpi=100)
	# plt.show()

def plot_histogram() : 
	""" plot histogram of number of nodes in each cell"""

	PATH = "/Users/jgzhang/Work/Densmore_lab/Partition/code_version/v2/genetic-circuit-partitioning/2021.4/runs/"
	sol_path = "results/electronic-circuits/alu/nparts/45/optimized_lc/"
	bm_file = PATH + "benchmark/electronic-circuits/alu/DAG.edgelist"

	G = load_graph (bm_file)
	in_nodes, out_nodes, nonprimitives  = gp.get_nonprimitive_nodes (G)
	G_primitive = gp.get_G_primitive (G, nonprimitives)

	solDict = gp.load_opt_part_sol (PATH + sol_path + 'part_solns.txt')
	iteration = 5
	part = solDict[iteration]['part']
	part_opt = [gp.get_part(part, n) for n in G_primitive.nodes()]
	# gp.visualize_assignment_graphviz (G, part_opt, nonprimitives, 'TRUE', PATH + sol_path, iteration, [])

	node_count = [len(part[p]) for p in part]
	print(node_count)
	x = sorted(set(node_count))
	y = [node_count.count(n) for n in x]

	fig = plt.figure(figsize=(5,5))
	ax = fig.add_subplot(111)
	ax.bar(x,y)
	ax.set_xlabel('No. nodes/cell')
	ax.set_ylabel('Count')
	plt.savefig(PATH+sol_path+str(iteration)+'_dist_nodecount.pdf', dpi=100)
	plt.show()



if __name__ == '__main__':
	PATH = '/home/ubuntu/genetic-circuit-partitioning/2021.4/'
	# PATH = '/Users/jgzhang/Work/Densmore_lab/Partition/code_version/v2/genetic-circuit-partitioning/2021.4/'
	

	# count_nodes (PATH + 'runs/results/4-input-boolean-circuits', PATH + 'runs/results/5-input-boolean-circuits')
	# partition_stats ()
	# compile_best_solutions ()
	# plot_timestep ()
	# count_motif()
	# plot_motif_occurence (PATH + 'runs/results/analysis/motif freq.txt')
	# plot_deltaD ()
	# G = load_graph (PATH + 'runs/benchmark/6-input-boolean-circuits/30/DAG.edgelist')
	# plot_network (G, PATH + 'runs/benchmark/bionetwork/random_graph/n50_p0.02/DAG.pdf')
	# plot_centrality (G, PATH + 'runs/benchmark/6-input-boolean-circuits/centrality.pdf')

	# best_soln = PATH + 'runs/results/5-input-boolean-circuits/best_solns.txt'
	# data = load_sol	(best_soln)
	# avgNode = np.mean([(int(data[bm]['Nodes'])-6)/int(data[bm]['Npart']) for bm in data])
	# print(avgNode)
	# partition = generate_edgelist_of_cells (PATH)
	# converted_names = convert_name (PATH + 'runs/results/electronic-circuits/md5Core/Jai_solution/gate_name_conversion.txt')
	visualize_subnetworks_unmet_constraint (PATH, 'lc')
	# compare_runs (PATH, 'lc')
	# compare_gatenum_distribution (PATH, 'lc')
	# compare_qs_distribution (PATH, 'lc')

	# plot_cell_edgelist (PATH+'runs/benchmark/electronic-circuits/alu/DAG.edgelist', PATH+'runs/results/electronic-circuits/alu/nparts/45/optimized_lc/part_iteration2.edgelist', partition)
	# g = nx.read_edgelist(PATH+'runs/benchmark/electronic-circuits/md5Core/DAG.edgelist', create_using=nx.DiGraph())
	# in_nodes, out_nodes, nonprimitives  = gp.get_nonprimitive_nodes (g)

	# parameter_scan_bionetworks (PATH)

	# plot_histogram ()
	
