""" this script plots the results of 4-input or 5-input boolean circuits """

import networkx as nx
import matplotlib.pyplot as plt 
import os
import numpy as np
from networkx.drawing.nx_agraph import graphviz_layout
from matplotlib.ticker import MaxNLocator
import pandas as pd
import copy
import itertools
from matplotlib import cm


def load_graph (dag_file):
	G = nx.read_edgelist(dag_file, nodetype = str, create_using=nx.DiGraph())
	return G

def load_part_sol (inputfile):
	lines = [open(inputfile, 'r').read().strip("\n")][0].split('\n')
	T = int( lines[0].split('\t')[1] )
	cut = int( lines[1].split('\t')[1] )
	partDict = {}
	for line in lines[2:]: 
		tokens = line.split('\t')
		part = int( tokens[0].split(' ')[-1] )
		nodes = tokens[1].split(',')
		partDict[part] = nodes
		for node in nodes: 
			if node in ['a', 'b', 'c', 'd', 'e']:
				print('partition error')
	return T, cut, partDict

def load_data (inputfile):
	lines = [open(inputfile, 'r').read().strip("\n")][0].split('\n')
	header = lines[0].split('\t')
	data_dict = {}
	for line in lines[1:]:
		tokens = line.split('\t')
		data_dict[tokens[0]] = {}
		for idx, token in enumerate(tokens, 1):
			data_dict[tokens[0]][header[idx-1]] = token
	return data_dict

def load_motif_data (inputfile):
 	lines = [open(inputfile, 'r').read().strip("\n")][0].split('\n')
 	motif_dict = {}
 	for line in lines[1:]:
 		tokens = line.split('\t')
 		motif_dict[tokens[0]] = float(tokens[1])
 	return motif_dict

def check_cycles (partG):
	""" check if partitioned cells contain cycles/loops """
	try:
		cycles = nx.find_cycle(partG)
		cycle_free = False
	except: 
		cycle_free = True
		pass 
	return cycle_free


def check_motif_allowed(matrix, stringent):
	""" check if all cells have allowed communications (2 in, 2 out at most)""" 
	motif_allowed = True
	out_deg = np.sum(matrix, axis=1)
	in_deg = np.sum(matrix, axis=0)
	summed_deg = np.sum(matrix, axis=1)+np.sum(matrix, axis=0)

	if stringent == True:
		if max(in_deg) > 2: # sum by rows 
			motif_allowed = False
		if max(out_deg) > 2: # sum by cols
			motif_allowed = False
		if max(summed_deg) > 4:
			motif_allowed = False
	else: 
		if max(summed_deg) > 4:
			motif_allowed = False
	return motif_allowed

def partition_matrix (DiG, partition, nonprimitives):

	# generate adjency matrix
	numParts = max(partition[1])+1
	matrix = np.zeros(shape=(numParts, numParts))

	G_primitive = copy.deepcopy(DiG)
	for node in nonprimitives:
		if node in G_primitive.nodes():
			G_primitive.remove_node(node)

	for edge in G_primitive.edges():
		v1, v2 = edge[0], edge[1]
		part_v1 = partition[1][list(G_primitive.nodes()).index(v1)]
		part_v2 = partition[1][list(G_primitive.nodes()).index(v2)]
		if part_v1 != part_v2:
			matrix[part_v1][part_v2] += 1

	# generate DAG representing cell-cell communication from the adjency matrix
	rows, cols = np.where(matrix != 0)
	edges = zip(rows.tolist(), cols.tolist())
	partG = nx.DiGraph() 
	partG.add_edges_from(edges) 

	return matrix, partG

def count_nodes (indir):
	""" count the number of nodes in graphs """
	count = []
	for i in range(1,100000):
		DiG = nx.DiGraph()
		dag_file = path+str(i)+'/DAG.edgelist'
		if os.path.exists(dag_file):
			DiG = nx.read_edgelist(dag_file)
			count.append( len(list(DiG.nodes())) )

	fig = plt.figure(figsize=(5,5)) 
	ax  = fig.add_subplot(111)
	ax.hist(count, ec='k', fc='wheat')
	ax.set_xlabel('Number of vertices')
	ax.set_ylabel('Count')
	# ax.set_xlim([0, 25])
	for spine in ['top', 'right']:
		ax.spines[spine].set_linewidth(0)
	plt.savefig('Node distribution - 5 input.pdf', dpi=200)
	plt.show()

def get_part (partDict, node): 
	for part, nodelist in partDict.items(): 
		if node in nodelist: 
			return part

def count_cuts (indir):
	""" count the number of cuts after partitioning circuits using METIS """
	cuts = []
	for i in range(1, 9001):
		pfile = path+str(i)+'/part_solns.txt'
		if os.path.exists(pfile):
			lines = [open(pfile, 'r').read().strip("\n")][0].split('\n')
			cut = lines[1].split('\t')[1]
			cuts.append(int(cut))

	x = np.arange(0, max(cuts)+1)
	occur = [cuts.count(c) for c in x]

	f_out = open('part occur - 5 input.txt', 'w')
	f_out.write('Cut\tOccur\n')
	for idx, c in enumerate(x): 
		f_out.write(str(c)+'\t'+str(occur[idx])+'\n')


def count_T (indir):
	""" count the number of cuts after partitioning circuits using METIS """
	cuts = []
	for i in range(1, 9001):
		pfile = path+str(i)+'/part_solns.txt'
		if os.path.exists(pfile):
			lines = [open(pfile, 'r').read().strip("\n")][0].split('\n')
			cut = lines[0].split('\t')[1]
			cuts.append(int(cut))

	x = np.arange(0, max(cuts)+1)
	occur = [cuts.count(c) for c in x]

	f_out = open('T occur - 5 input.txt', 'w')
	f_out.write('T\tOccur\n')
	for idx, c in enumerate(x): 
		f_out.write(str(c)+'\t'+str(occur[idx])+'\n')


def count_opt_params (indir):
	""" count the number of T and cut after METIS and module shifting"""
	cuts = []
	cuts_o = []
	Ts = []
	Ts_o = []
	cut_opt_num = 0
	T_opt_num = 0
	cut_T_opt_num = 0
	f_out = open('./boolean circuit/optimization improvement - 4 input (112220).txt', 'w')
	f_out.write('Graph\tNodes\tValid Motif_METIS_lowconstraint\tValid Motif_METIS_highconstraint\tValid Motif_Optimized_lowconstraint\tValid Motif_Optimized_highconstraint\tCycle Free_METIS\tCycle Free_Optimized\tT_Metis\tT_Optimized\tT_Delta\tCut_Metis\tCut_Optimized\tCut_Delta\n')

	nonprimitives = ['a', 'b', 'c', 'd', 'e', 'out']

	for i in range(1, 3000):
		print(i)

		G = load_graph (path+str(i)+'/DAG.edgelist')
		G_primitive = copy.deepcopy(G)
		for node in nonprimitives:
			if node in G.nodes():
				G_primitive.remove_node(node)

		T_o, cut_o, partDict_o = load_part_sol (path+str(i)+'/part_solns.txt')

		# check if original partition modules satisfy constraints
		part_opt_format = (cut_o, [get_part(partDict_o, n) for n in G_primitive.nodes()])
		matrix, partG = partition_matrix (G, part_opt_format, nonprimitives)
		cycle_free_o = check_cycles(partG)
		motif_valid_o_lc, motif_valid_o_hc = check_motif_allowed(matrix, False), check_motif_allowed(matrix, True)

		Ts_o.append(T_o)
		cuts_o.append(cut_o)
		node_num = len (sum([partDict_o[part] for part in partDict_o], []))
		cut_opt = False
		T_opt = False
		cut_T_opt = False
		minT, minC = T_o, cut_o
		motif_valid_lc, motif_valid_hc, cycle_free = 'NA', 'NA', 'NA'
		
		# if graph i has an optimized solution 
		if os.path.exists(path+str(i)+'/optimized_hc/part_opt1.txt'):
			cycle_free, motif_valid_hc = True, True
			optfiles = [filename for filename in os.listdir(path+str(i)+'/optimized_hc') if filename.startswith('part_opt') and filename.endswith('.txt')]
			for optfile in optfiles:
				T, cut, partDict = load_part_sol (path+str(i)+'/optimized_hc/'+optfile)
				part_opt_format = (cut, [get_part(partDict, n) for n in G_primitive.nodes()])
				matrix, partG = partition_matrix (G, part_opt_format, nonprimitives)
				if T < minT: minT = T
				if cut < minC: minC = cut
		elif os.path.exists(path+str(i)+'/optimized_lc/part_opt1.txt'):
			print('low constraint')
			cycle_free, motif_valid_lc = True, True
			optfiles = [filename for filename in os.listdir(path+str(i)+'/optimized_lc') if filename.startswith('part_opt') and filename.endswith('.txt')]
			for optfile in optfiles:
				T, cut, partDict = load_part_sol (path+str(i)+'/optimized_lc/'+optfile)
				part_opt_format = (cut, [get_part(partDict, n) for n in G_primitive.nodes()])
				matrix, partG = partition_matrix (G, part_opt_format, nonprimitives)
				if T < minT: minT = T
				if cut < minC: minC = cut

		cuts.append(minC)
		Ts.append(minT)

		f_out.write('\t'.join([str(i), str(node_num), str(motif_valid_o_lc), str(motif_valid_o_hc), str(motif_valid_lc), str(motif_valid_hc), str(cycle_free_o), str(cycle_free), str(T_o), str(minT), str(T_o-minT), str(cut_o), str(minC), str(cut_o-minC)])+'\n')
	f_out.close()

	# write occurrence of cuts and T before optimization 
	x = np.arange(0, max(cuts_o)+1)
	occur = [cuts_o.count(c) for c in x]

	f_out = open('./boolean circuit/cut occur before optimization - 4 input (G1-3000).txt', 'w')
	f_out.write('Cut\tOccur\n')
	for idx, c in enumerate(x): 
		f_out.write(str(c)+'\t'+str(occur[idx])+'\n')
	f_out.close()

	x = np.arange(0, max(Ts_o)+1)
	occur = [Ts_o.count(t) for t in x]

	f_out = open('./boolean circuit/T occur before optimization - 4 input (G1-3000).txt', 'w')
	f_out.write('T\tOccur\n')
	for idx, t in enumerate(x): 
		f_out.write(str(t)+'\t'+str(occur[idx])+'\n')
	f_out.close()

	# write occurence of cuts and T after optimization 
	x = np.arange(0, max(cuts)+1)
	occur = [cuts.count(c) for c in x]

	f_out = open('./boolean circuit/cut occur after optimization - 4 input (G1-3000).txt', 'w')
	f_out.write('Cut\tOccur\n')
	for idx, c in enumerate(x): 
		f_out.write(str(c)+'\t'+str(occur[idx])+'\n')
	f_out.close()

	x = np.arange(0, max(Ts)+1)
	occur = [Ts.count(t) for t in x]

	f_out = open('./boolean circuit/T occur after optimization - 4 input (G1-3000).txt', 'w')
	f_out.write('T\tOccur\n')
	for idx, t in enumerate(x): 
		f_out.write(str(t)+'\t'+str(occur[idx])+'\n')	
	f_out.close()


def count_cell_module (part_matrix_file):
	moduleList = []
	matrix = np.load(part_matrix_file)
	# print(matrix)
	summed_row = np.sum(matrix, axis=1).tolist() # sum by each row
	summed_col = np.sum(matrix, axis=0).tolist() # sum by each col
	# print(summed_row)
	# print(summed_col)
	for cell in range(matrix.shape[0]):
		module = (summed_col[cell], summed_row[cell]) # (in degree, out degree)
		moduleList.append(module)
	return moduleList



def get_cell_comm_module_opt (path): 
	moduleDict = {}
	for i in range(1, 65530):
		print(i)
		file_path = path + str(i)
		dag_file = path + '/DAG.edgelist'
		G = load_graph (dag_file)
		# part_sol = path + '/part_solns.txt'
		mfile =  file_path + '/part_matrix.npy'

		if os.path.exists(file_path + '/optimized_hc/part_opt1.npy'):
			moduleList = count_cell_module(file_path + '/optimized_hc/part_opt1.npy')
		elif os.path.exists(file_path + '/optimized_lc/part_opt1.npy'):
			moduleList = count_cell_module(file_path + '/optimized_lc/part_opt1.npy')
		else: 
			moduleList = count_cell_module(mfile)
		print(moduleList)
		# moduleList = count_cell_module(G, mfile)

		for module in moduleList:
			if module in moduleDict:
				moduleDict[module] += 1
			else: 
				moduleDict[module] = 1

	# transfer module number to a matrix
	dim = 0
	for key in moduleDict.keys():
		indegree, outdegree = key[0], key[1]
		if max([indegree, outdegree]) > dim:
			dim = max([indegree, outdegree])
	modulematrix = np.zeros(shape=(int(dim+1), int(dim+1)))
	for indegree in range(int(dim+1)):
		for outdegree in range(int(dim+1)):
			module = (indegree, outdegree)
			if module in moduleDict:
				modulematrix[indegree][outdegree] = moduleDict[module]
	print(modulematrix)
	# np.save('./boolean circuit/module matrix METIS (G1-65529).npy', modulematrix)



def plot_cell_comm_module (infile):

	modulematrix = np.load(infile)
	dim = modulematrix.shape[0]
	print(pd.DataFrame(modulematrix))

	# plot module x - num of in edges, y - num of out edges
	fig = plt.figure(figsize=(5,5))
	ax = fig.add_subplot(111)
	df = pd.DataFrame(modulematrix)
	# log transform count
	df = df.apply(lambda x: np.log(x))
	print(df)
	# df[5] = [0,0,0,0,0]   # append a col
	# df.loc[5] = [0,0,0,0,0,0].  # append a row 
	# print(df)
	# scale count from 0 to 1
	# df = (df - min(moduleList.values())) /(max(moduleList.values()) - min(moduleList.values()))

	plt.imshow(df, cmap='Greys', vmin=0, vmax=7.5)
	plt.colorbar()
	ax.set_yticks(list(range(int(dim))))
	ax.set_xticks(list(range(int(dim))))
	# ax.set_xticks(list(range(6)))
	# ax.set_yticks(list(range(6)))
	
	plt.savefig('./boolean circuit/module frequency METIS - ln (G1-3000).pdf', dpi=200)
	plt.show()


def to_tuple (matrix):
	return tuple(tuple(arr.tolist()) for arr in matrix)


def get_motif_occurence (path):
	""" plot different cell-cell communication motif """
	motifDict = {}
	for i in range(1, 65530):
		print('graph', i)
		file_path = path + str(i)
		dag_file = path + '/DAG.edgelist'
		G = load_graph (dag_file)
		# part_sol = path + '/part_solns.txt'
		mfile =  file_path + '/part_matrix.npy'

		if os.path.exists(mfile):

			if os.path.exists(file_path + '/optimized_hc/part_opt1.npy'):
				optfiles = [path+str(i)+'/optimized_hc/'+filename for filename in os.listdir(path+str(i)+'/optimized_hc') if filename.startswith('part_opt') and filename.endswith('.npy')]
			elif os.path.exists(file_path + '/optimized_lc/part_opt1.npy'):
				optfiles = [path+str(i)+'/optimized_lc/'+filename for filename in os.listdir(path+str(i)+'/optimized_lc') if filename.startswith('part_opt') and filename.endswith('.npy')]
			else: 
				optfiles = [file_path+'/part_matrix.npy']
			# print(optfiles)

			for optfile in optfiles:
				matrix = np.load(optfile)
				dim = matrix.shape[0]
				# print('part matrix', matrix)
				I = np.identity(n=dim)  # identity matrix
				dmList = []
				for pm in itertools.permutations(I):  # permutation matrix
					pm = np.array(pm)
					# print('perm matrix', pm)
					dm = np.dot(pm, matrix)               # permutate by rows
					# print('product', np.dot(dm, pm.T))
					dm = np.dot(dm, pm.T)
					dmList.append (to_tuple(dm))              # append all diagnizable matrix of partition matrix 
				# print('all diag matrix', dmList)

				# motif_exist = set(dmList).intersect(set(motifDict.keys()))
				if any(dm in motifDict for dm in dmList) == False:
					# print('motif not found in dict', matrix)
					motifDict[to_tuple(matrix)] = 1/len(optfiles)
					# print('new motifdict', motifDict)
				else:
					motif = set(dmList).intersection(set(motifDict.keys()))
					# print('motif', motif,'already in dict')
					motifDict[tuple(motif)[0]] += 1/len(optfiles)


	# write output 
	f_out = open('./boolean circuit/motif freq.txt', 'w')
	f_out.write('Motif matrix\tOccur\n')
	for motif in motifDict:
		f_out.write(str(motif)+'\t'+str(motifDict[motif])+'\n')
		print(motif, motifDict[motif])

def to_array (strel):
	return np.array(eval(strel))

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
	fig = plt.figure(figsize=(12,12))
	sp = 1
	for motif in sort_occur:
		print(sp, motif)
		ax = fig.add_subplot(12,12,sp)
		matrix = np.array(eval(motif))
		partG = generate_comm_graph(matrix)
		pos = graphviz_layout(partG, prog='dot')
		nx.draw_networkx_nodes(partG, pos, node_size=30, node_color='#b18ea6')
		labels={n:n for n in partG.nodes()}
		nx.draw_networkx_edges(partG, pos)
		ax.axis('off')
		sp += 1

	plt.savefig('./boolean circuit/motif occur - test.pdf', dpi=200)
	plt.show()

	## plot bar plot of motif occur
	# fig = plt.figure(figsize=(5,3))
	# ax = fig.add_subplot(111)
	# x = np.arange(1, len(sort_occur_freq)+1)
	# labels = ['M'+str(i) for i in range(1, len(x)+1)]
	# sizes = [100*i for i in sort_occur_freq]
	# explode = (0.1,)*len(x)
	# ax.pie(sizes, explode=explode, labels=labels, autopct='%1.1f%%',
 #        shadow=False, startangle=0)
	# plt.savefig('./boolean circuit/motif occur freq - pie.pdf', dpi=200)
	# plt.show()


def calc_opt_results (inputfile):
	""" after optimization, calculate how many graphs have solutions """
	data = load_data (inputfile)

	motif_valid_ori_lc = [data[graph]['Valid Motif_METIS_lowconstraint'] for graph in data].count('True')* 100 / len(data)
	motif_valid_ori_hc = [data[graph]['Valid Motif_METIS_highconstraint'] for graph in data].count('True')* 100 / len(data)
	motif_valid_opt_lc = [data[graph]['Valid Motif_Optimized_lowconstraint'] for graph in data].count('True')* 100 / len(data)
	motif_valid_opt_hc = [data[graph]['Valid Motif_Optimized_highconstraint'] for graph in data].count('True') * 100 / len(data)
	cycle_free_ori = [data[graph]['Cycle Free_METIS'] for graph in data].count('True') * 100 / len(data)
	cycle_free_opt = [data[graph]['Cycle Free_Optimized'] for graph in data].count('True')* 100 / len(data)

	# Solutions with valid modules after node migration
	sols_module_valid = 0
	for graph in data:
		module_valid = data[graph]['Valid Motif_METIS_lowconstraint'] 
		if module_valid == 'False':
			if data[graph]['Valid Motif_Optimized_lowconstraint'] == 'True':
				sols_module_valid += 1
		else: 
			sols_module_valid += 1
	print(100*sols_module_valid/len(data))

	sols_cycle_valid = 0
	for graph in data:
		cycle_valid = data[graph]['Cycle Free_METIS']
		if cycle_valid == 'False':
			if data[graph]['Cycle Free_Optimized'] == 'True':
				sols_cycle_valid += 1
		else: 
			sols_cycle_valid += 1
	print(100*sols_cycle_valid/len(data))

	print(motif_valid_ori_lc, motif_valid_ori_hc, motif_valid_opt_lc, motif_valid_opt_hc, cycle_free_ori, cycle_free_opt)


def compare_cuts (infile1, infile2): 
	""" compare cuts before and after optimization """
	lines1 = [open(infile1, 'r').read().strip("\n")][0].split('\n')
	lines2 = [open(infile2, 'r').read().strip("\n")][0].split('\n')
	cutDict = {}
	for line in lines1[1:]:
		tokens = line.split('\t')
		cutDict[int(tokens[0])] = {'metis': int(tokens[1])}

	for line in lines2[1:]:
		tokens = line.split('\t')
		cutDict[int(tokens[0])]['optimized'] = int(tokens[1])

	print(cutDict)
	x1, x2 = [], []
	y1, y2 = [], []
	for i in range(max(cutDict.keys())+1):
		x1.append( i-0.1 )
		x2.append( i+0.1 )
		y1.append(cutDict[i]['metis'])
		try:
			y2.append(cutDict[i]['optimized'])
		except KeyError:
			y2.append(0)

	### plot comparison 
	fig = plt.figure(figsize=(5,5))
	ax = fig.add_subplot(111)
	# ax.set_xlim([0, 7])
	ax.xaxis.set_major_locator(MaxNLocator(integer=True))
	(markerline, stemlines, baseline) = plt.stem(x1, y1, linefmt='#94b4a4', label='Initial Graph Partition')
	(markerline2, stemlines2, baseline2) = plt.stem(x2, y2, linefmt='#c56183', label='Node Migration Optimized')
	ax.set_xlabel('Cuts')
	ax.set_ylabel('Count')
	for spine in ['top', 'right']:
		ax.spines[spine].set_linewidth(0)
	plt.setp(baseline, visible=False)
	plt.setp(baseline2, visible=False)
	# plt.legend()
	markerline.set_markerfacecolor('#94b4a4')
	markerline2.set_markerfacecolor('#c56183')
	plt.savefig('./boolean circuit/T count comparison - 4 input (G1-3000).pdf', dpi=200)
	plt.show()	


def get_moved_node_property (G, dict1, dict2):
	""" compare partition dict to 
	1.  calculate how many nodes are moved to a different cell 
	2.  averaged moved distance """
	num_moved_nodes = 0
	for node in G.nodes():
		part_o = get_part (dict1, node)
		part_n = get_part (dict2, node)
		# eccentricity = nx.eccentricity(G, node)
		
		if part_o != part_n: 
			num_moved_nodes += 1

	return num_moved_nodes 

def calc_distance_to_boundry (G, part_dict, node, p):
	""" for a given node, calculate its closest distance to the boundry nodes within the block """
	block = get_part(part_dict, node) # the block where this node is located 
	# find the nodes at the boundry of this block (that has connection to nodes in other blocks)
	distanceList = []   # distance to nodes at the boundry
	for n in nodes: 
		nbs = list(G.neighbors(n))      # neighboring nodes of each node in this block
		others = set(nbs) - set(nodes)  # others: neighboring nodes that are not in the same block
		if len(others) != 0:
			dist = len(p[n][node])-1    
			distanceList.append(dist)
	# return minimal distance to the boundry 
	return min(distanceList)   

def plot_distance_to_boundry (data):
	""" plot the distance to block boundary of moved nodes """
	dist_dict, opt_dist_dict, opt_dist_dict2 = {}, {}, {}
	for n in data:
		dist_dict[n] = {}
		opt_dist_dict[n] = {}
		for g in data[n]:
			G = nx.read_edgelist('./Random graphs/n'+str(n)+'/G'+str(n)+'.'+str(g)+'.edgelist', nodetype = int)
			p = nx.shortest_path(G)  # shortest path connecting two nodes in this graph 
			
			if min(data[n][g]) not in [0, -1]: 
				mindelta = min([data[n][g][s]['New Cut']-data[n][g][s]['Ori Cut'] for s in data[n][g]])
				for s in data[n][g]:
					dcut = data[n][g][s]['New Cut']-data[n][g][s]['Ori Cut']
					part_ori = data[n][g][s]['Ori Partition Result']
					part_new = data[n][g][s]['New Partition']

					for idx in range(len(part_ori)):
						if part_ori[idx] != part_new[idx]:
							mindist = calc_distance_to_boundry (G, part_ori, idx, p)
							if mindist not in dist_dict[n]:
								dist_dict[n][mindist] = 1
							else: 
								dist_dict[n][mindist] += 1

					if dcut == mindelta:
						for idx in range(len(part_ori)):
							if part_ori[idx] != part_new[idx]:
								mindist = calc_distance_to_boundry (G, part_ori, idx, p)
								if mindist not in opt_dist_dict[n]:
									opt_dist_dict[n][mindist] = 1
								else: 
									opt_dist_dict[n][mindist] += 1

	print(dist_dict)
	print(opt_dist_dict)
	dist_dict = convert_count_to_perc(dist_dict)
	opt_dist_dict = convert_count_to_perc(opt_dist_dict)
	print(dist_dict)
	print(opt_dist_dict)

	fig = plt.figure(figsize=(8,5))
	ax = fig.add_subplot(121)
	for n in range(20, 60, 10): 
		plist = dist_dict[n].values()
		ax.plot(plist, 'o-', label=n)
	#ax.legend()
	ax.set_ylabel('Frequency', fontsize=12)
	ax.set_ylim([0,100])
	for spine in ['top', 'right']:
		ax.spines[spine].set_linewidth(0)

	ax = fig.add_subplot(122)
	for n in range(20, 60, 10): 
		plist = opt_dist_dict[n].values()
		ax.plot(plist, 'o-', label=n)
	ax.set_ylim([0,100])
	ax.legend()

	plt.xlabel('Closest distance to boundry', fontsize=12)
	for spine in ['top', 'right']:
		ax.spines[spine].set_linewidth(0)
	plt.savefig('./data/distance to boundry.pdf')
	plt.show()


def count_num_moved_nodes (path):
	""" in each optimized graph, count how many nodes are moved to a different cell """

	nonprimitives = ['a', 'b', 'c', 'd', 'e', 'out']
	num_moved_nodes_dict = {}

	f_out = open ('./boolean circuit/number of moved nodes.txt', 'w')
	f_out.write('Graph\tNumber of moved nodes in solution\n')

	for i in range(1, 65530):

		if os.path.exists(path+str(i)+'/DAG.edgelist'):
			G = load_graph (path+str(i)+'/DAG.edgelist')
			G_primitive = copy.deepcopy(G)
			for node in nonprimitives:
				if node in G.nodes():
					G_primitive.remove_node(node)
			
			try:
				T_o, cut_o, partDict_o = load_part_sol (path+str(i)+'/part_solns.txt')

				# check if original partition modules satisfy constraints
				# part_opt_format = (cut_o, [get_part(partDict_o, n) for n in G_primitive.nodes()])
				# matrix, partG = partition_matrix (G, part_opt_format, nonprimitives)

				# if graph i has an optimized solution 
				if os.path.exists(path+str(i)+'/optimized_hc/part_opt1.txt'):
					# print('high constraint')
					optfiles = [filename for filename in os.listdir(path+str(i)+'/optimized_hc') if filename.startswith('part_opt') and filename.endswith('.txt')]
					num_moved_nodes_i = []
					for optfile in optfiles:
						T, cut, partDict = load_part_sol (path+str(i)+'/optimized_hc/'+optfile)
						# part_opt_format = (cut, [get_part(partDict, n) for n in G_primitive.nodes()])
						# matrix, partG = partition_matrix (G, part_opt_format, nonprimitives)
						num_moved_nodes = get_moved_node_property (G_primitive, partDict_o, partDict)
						num_moved_nodes_i.append (str(num_moved_nodes))
					
					f_out.write(str(i)+'\t'+str(len(partDict.keys()))+'\t'+','.join(num_moved_nodes_i)+'\n')

				elif os.path.exists(path+str(i)+'/optimized_lc/part_opt1.txt'):
					# print('low constraint')
					optfiles = [filename for filename in os.listdir(path+str(i)+'/optimized_lc') if filename.startswith('part_opt') and filename.endswith('.txt')]
					num_moved_nodes_i = []
					for optfile in optfiles:
						T, cut, partDict = load_part_sol (path+str(i)+'/optimized_lc/'+optfile)
						# part_opt_format = (cut, [get_part(partDict, n) for n in G_primitive.nodes()])
						# matrix, partG = partition_matrix (G, part_opt_format, nonprimitives)
						num_moved_nodes = get_moved_node_property (G_primitive, partDict_o, partDict)
						num_moved_nodes_i.append (str(num_moved_nodes))

					f_out.write(str(i)+'\t'+str(len(partDict.keys()))+'\t'+','.join(num_moved_nodes_i)+'\n')
			except IndexError:
				pass 



def plot_num_moved_nodes (inputfile):

	num_moved_nodes_dict = dict()
	lines = [open(inputfile, 'r').read().strip("\n")][0].split('\n')
	for line in lines[1:]:
		num = [int(n) for n in (line.split('\t')[2]).split(',')]
		part  = int(line.split('\t')[1])
		if part not in num_moved_nodes_dict:
			num_moved_nodes_dict[part] = {}
		for n in num:
			if n not in num_moved_nodes_dict[part]:
				num_moved_nodes_dict[part][n] = 1
			else: 
				num_moved_nodes_dict[part][n] += 1
	print(num_moved_nodes_dict)
	fig = plt.figure(figsize=(8,5))
	ax  = fig.add_subplot(111)
	x = np.arange(5)
	y1 = [v/sum([num_moved_nodes_dict[2][n] for n in range(1, 6)]) for v in [num_moved_nodes_dict[2][n] for n in range(1, 6)]]
	y2 = [v/sum([num_moved_nodes_dict[3][n] for n in range(1, 6)]) for v in [num_moved_nodes_dict[3][n] for n in range(1, 6)]]
	y3 = [v/sum([num_moved_nodes_dict[4][n] for n in range(1, 6)]) for v in [num_moved_nodes_dict[4][n] for n in range(1, 6)]]
	width = 0.2

	# plot data in groups
	plt.bar(x-0.2, y1, width, color='#ecb390')
	plt.bar(x, y2, width, color='#ecdfc8')
	plt.bar(x+0.2, y3, width, color='#df7861')
	plt.legend(['2 Cells', '3 Cells', '4 Cells'])
	plt.xticks(x, [1,2,3,4,5])
	plt.ylim([0, 0.5])
	for spine in ['top', 'right']:
		ax.spines[spine].set_linewidth(0)
	ax.set_xlabel('Number of moved nodes')
	ax.set_ylabel('Count')
	plt.subplots_adjust(left=0.17, bottom=0.15)
	plt.savefig('Number of moved nodes.pdf', dpi=200)
	plt.show()






if __name__ =='__main__':
	# path = '/Users/jgzhang/Programs/Cello2/sample-input/DNACompiler/4-input/'
	# load_data ('/Users/jgzhang/Work/Densmore lab/Partition/boolean circuit/optimization improvement - 4 input (111820).txt')
	# count_nodes(path)
	# plot_modules(path)
	# count_opt_cuts(path)
	# cal_delta_cuts (path)
	# plot_delta_cut ('./boolean circuit/optimization improvement - 4 input.txt')
	# plot_cuts('./boolean circuit/part occur after optimization - 4 input.txt')
	# compare_cuts ('./boolean circuit/T occur before optimization - 4 input (G1-3000).txt', './boolean circuit/T occur after optimization - 4 input (G1-3000).txt')
	# get_cell_comm_module_opt (path)
	# plot_cell_comm_module ('./boolean circuit/module matrix METIS (G1-3000).npy')
	# count_cycles(path)
	# count_opt_params (path)
	# calc_opt_results('/Users/jgzhang/Work/Densmore lab/Partition/boolean circuit/optimization improvement - 4 input (112220).txt')
	# get_motif_occurence (path)
	# plot_motif_occurence ('./boolean circuit/motif freq.txt')
	# count_num_moved_nodes (path)
	# plot_num_moved_nodes ('./boolean circuit/number of moved nodes.txt')


