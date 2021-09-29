import networkx as nx
import numpy as np
import random
import copy
import math
import os
import itertools
import metis
from networkx.drawing.nx_agraph import graphviz_layout
import matplotlib.pyplot as plt
from matplotlib import gridspec
import shutil
import seaborn as sns
import warnings
warnings.filterwarnings("ignore")

def load_file (inputfile):
	lines = [open(inputfile, 'r').read().strip("\n")][0].split('\n')
	header = lines[0].split('\t')
	data = {}
	for line in lines[1:]:
		tokens = line.split('\t')
		data[tokens[0]] = {}
		for idx, token in enumerate(tokens[1:]):
			data[tokens[0]][header[idx+1]] = token
	return data


def load_graph (dag_file):
	G = nx.read_edgelist(dag_file, nodetype = str, create_using=nx.DiGraph())
	innodes, outnodes = [], []
	for node in G.nodes():
		indegree = G.in_degree(node)
		outdegree = G.out_degree(node)
		if outdegree == 0:
			outnodes.append(node)
		if indegree == 0:
			innodes.append(node)
	return G, innodes, outnodes


def load_ori_part_sol (inputfile):
	lines = [open(inputfile, 'r').read().strip("\n")][0].split('\n')
	cut = int( lines[0].split('\t')[1] )
	partDict = {}
	for line in lines[1:]: 
		tokens = line.split('\t')
		part = int( tokens[0].split(' ')[-1] )
		nodes = tokens[1].split(',')
		partDict[part] = nodes
		# for node in nodes: 
		# 	if node in ['a', 'b', 'c', 'd', 'e']:
		# 		print('partition error')
	return cut, partDict

def load_part_sol (inputfile):
	lines = [open(inputfile, 'r').read().strip("\n")][0].split('\n')
	cut = int( lines[0].split('\t')[1] )
	partDict = {}
	for line in lines[1:]: 
		tokens = line.split('\t')
		part = int( tokens[0].split(' ')[-1] )
		nodes = tokens[1].split(',')
		partDict[part] = nodes
		# for node in nodes: 
		# 	if node in ['a', 'b', 'c', 'd', 'e']:
		# 		print('partition error')
	return cut, partDict



def load_opt_part_sol (inputfile):
	lines = [open(inputfile, 'r').read().strip("\n")][0].split('\n')
	iteration_idx = [lines.index(line) for line in lines if line.startswith('path')]
	print('interation index', iteration_idx)
	partDict, TDict, cutDict = {}, {}, {}
	for i, idx in enumerate(iteration_idx):
		print('i-th', i)
		print('index', idx)
		iteration = int (lines[idx].split('\t')[1] )
		T = int( lines[idx+1].split('\t')[1] )
		cut = int( lines[idx+2].split('\t')[1] )
		partDict[iteration] = {}
		TDict[iteration] = T
		cutDict[iteration] = cut
		if idx!= iteration_idx[-1] and len(iteration_idx)!= 1: 
			part_lines = lines[idx+3:iteration_idx[i+1]]
			# print(part_lines)
		else: 
			part_lines = lines[idx+3:]
			# print(part_lines)
		for line in part_lines: 
			tokens = line.split('\t')
			part = int( tokens[0].split(' ')[-1] )
			nodes = tokens[1].split(',')
			partDict[iteration][part] = nodes
	print(partDict)
	return TDict, cutDict, partDict


def write_part_sol (path, partDict, T, cut):
	f_out = open(path+'/part_solns.txt', 'w')
	f_out.write('T\t'+str(max(T))+'\n')
	f_out.write('cut\t'+str(cut)+'\n')
	for part in range(max(partDict)+1):
		f_out.write('Partition '+str(part)+'\t')
		f_out.write(','.join([str(n) for n in partDict[part]])+'\n')


def load_optimized_sol (inputfile):
	lines = [open(inputfile, 'r').read().strip("\n")][0].split('\n')
	solDict = {}
	start_idxs = [n for n, l in enumerate(lines) if l.startswith('T')]
	for idx, start_idx in enumerate(start_idxs):
		solDict[idx] = {}
		if start_idx != start_idxs[-1]:
			end_idx = start_idxs[idx+1]
		else: 
			end_idx = len(lines)
		sol_lines = lines[start_idx:end_idx]
		T = sol_lines[0].split('\t')[1]
		cut = sol_lines[1].split('\t')[1]
		solDict[idx]['T'] = int(T)
		solDict[idx]['cut'] = int(cut)
		solDict[idx]['part'] = {}
		for l in sol_lines[2:]:
			print(l)
			if l.startswith('path') == False:
				part = int((l.split('\t')[0]).split(' ')[1])
				nodes = (l.split('\t')[1]).split(',')
				solDict[idx]['part'][part] = nodes 
	return solDict

def load_minT_sol (inputfile):
	lines = [open(inputfile, 'r').read().strip("\n")][0].split('\n')
	minTDict = {}
	for line in lines[1:]:
		iteration, minT = int(line.split('\t')[0]), line.split('\t')[1]
		minTDict[iteration] = minT 
	return minTDict



# def load_opt_part_sol (inputfile):

# 	lines = [open(inputfile, 'r').read().strip("\n")][0].split('\n')
# 		# if graph i has an optimized solution
# 		# pfile = path+str(i)+'/optimized_solution_'+str(i)+'.txt'
# 		# if os.path.exists(pfile):
# 	minCut = 100
# 	for line in lines:
# 		if line.startswith('cut'):
# 			cut = int( lines[0].split('\t')[1] )
# 			if cut < minCut:
# 				minCut = cut
			

def assign_qs (dag_file, part_matrix, part_sol):
	partG = load_part_graph (part_matrix)
	cut, partDict = load_part_sol (part_sol)
	# a dictionary recording qs used in a given partition 
	qsDict = {}
	for part in partDict: 
		qsDict[part] = []

	for edge in list(partG.edges()):
		print('edge', edge)
		p_s, p_e = edge[0], edge[1]  # partition where the edge begins and ends
		qs_used = qsDict[p_e] + qsDict[p_s]  # all the qs that have been used 
		print('used qs in', p_s, p_e, qs_used)
		p_o = list(partG.nodes())    # partitions other than p_s and p_e
		p_o.remove(p_s)
		p_o.remove(p_e)
		for p in p_o:
			paths = list( nx.all_simple_paths(partG, source=p_s, target=p) )
			print('all paths from', p_s, 'to p', p, paths)
			p_in_paths = set([part for path in paths for part in path])
			for p_p in p_in_paths:
				qs_used.append(qsDict[p_p])
				print(p_p, ' in path', 'has qs', qsDict[p_p])
		print('qs used by downstream nodes', set(qs_used))
		# assign qs to edge 
		if len(set(qs_used)) == 0:
			partG.edges[p_s, p_e, 0]['qs'] = 1
			print('assigning qs to edge', edge, 1)
		else: 
			partG.edges[p_s, p_e, 0]['qs'] = max(set(qs_used)) + 1
			print('assigning qs to edge', edge, max(set(qs_used))+1 )

		# update qsDict 
		qsDict[p_s].append(partG.edges[p_s, p_e]['qs'])
		qsDict[p_e].append(partG.edges[p_s, p_e]['qs'])



def load_part_graph (inputfile):
	G = nx.MultiDiGraph()
	matrix = np.load(inputfile)
	for i in range(matrix.shape[0]):
		for j in range(matrix.shape[1]):
			num_edge = matrix[i][j]  # number of edges from i to j
			for e in range(int(num_edge)):
				G.add_edge(i, j, key=0)
	print('edges', list(G.edges()))
	return G


def get_part (partDict, node): 
	for part, nodelist in partDict.items(): 
		if node in nodelist: 
			return part

 
def cal_cut (G, partDict):
	""" calculate the cut in a new partition """
	cut = 0
	for edge in G.edges():
		part0 = get_part(partDict, edge[0])
		part1 = get_part(partDict, edge[1])
		if part0 != part1: 
			cut += 1
	return cut




def classify_nodes_subG (G, nodes):
	""" categorize nodes in G to innodes, outnodes of subgraph of nodes"""
	innodes, outnodes = [], []
	for e in G.edges():
		from_e, to_e = e[0], e[1]
		if to_e in nodes and from_e not in nodes:
			innodes.append(from_e)
		if from_e in nodes and to_e not in nodes:
			outnodes.append(to_e)
	# if input nodes/output nodes are in the subgraph 
	for n in nodes:
		if G.in_degree(n) == 0:
			innodes.append(n)
	for n in nodes:
		if G.out_degree(n) == 0:
			outnodes.append(n)
	return innodes, outnodes


def generate_combinations (n, rlist):
	""" from n choose r elements """
	combs = [list(itertools.combinations(n, r)) for r in rlist]
	combs = [item for sublist in combs for item in sublist]
	return combs


def calc_signal_path (dag, innodes, outnodes, partDict):
	""" count the number of boundaries that each input has to pass in order to reach the output """
	crosslist = []
	for inode in innodes:
		for outnode in outnodes:
			for path in sorted(nx.all_simple_edge_paths(dag, inode, outnode)):
				if (all(e in list(dag.edges()) for e in path)): # valid path that makes a directed path from source to target
					cross = 0
					for e in path:
						if (any(n in innodes+outnodes for n in e)) == False:
							if get_part(partDict, e[0]) != get_part(partDict, e[1]):
								cross += 1
					crosslist.append(cross)
	return crosslist


def calc_signal_path_subgraph (dag, innodes, outnodes, subG_nodes, partDict):
	""" count the number of boundaries that each input has to pass in order to reach the output """
	crosslist = []
	for inode in innodes:
		for outnode in outnodes:
			for path in sorted(nx.all_simple_edge_paths(dag, inode, outnode)):
				if (all(e in list(dag.edges()) for e in path)): # valid path that makes a directed path from source to target
					cross = 0
					for e in path:
						if (any(n in subG_nodes for n in e)) == True:
							if get_part(partDict, e[0]) != get_part(partDict, e[1]):
								cross += 1
					crosslist.append(cross)
	return crosslist



def compare_crosses (dag, partDict, partDict2):
	cross = calc_signal_path (dag, partDict)
	cross_opt = calc_signal_path (dag, partDict2)
	print(cross)
	print(cross_opt)


def check_loops (matrix):
	""" check if partitioned cells form any loops """
	loop_free = True
	for i in range(matrix.shape[0]):
		for j in range(matrix.shape[1]):
			if i != j:
				if 0 not in [matrix[i][j], matrix[j][i]]:
					loop_free = False
	return loop_free

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

   
def optimize_signal(G, innodes, outnodes, cut, partDict, primitive_only, maxNodes, motif_constraint, outdir):
	""" 
	optimize based on signal travel time from inputs to the output 
	1. calculate the times that inputs have to traverse cell boundries 
	2. optmize the max traverse to be as low as possible 
	"""

	Smin, Smax = 3, 7

	nonprimitives = innodes + outnodes
	G_primitive = copy.deepcopy(G)
	if primitive_only:
		for node in nonprimitives:
			if node in G.nodes():
				G_primitive.remove_node(node)

	T = calc_signal_path (G, innodes, outnodes, partDict)
	minT = max(T)
	solN = 0
	print('original traverse time', T)

	### optimize traverse times 
	# make a directory to store results
	if os.path.exists(outdir):
		shutil.rmtree(outdir)
		os.mkdir(outdir)
	else: 
		os.mkdir(outdir)

	solfile = outdir+'sigaltraverse_optimized_solution.txt'
	potentialsolfile = outdir+'potential_solution.txt'

	if minT > 0:
		# choose nodes to move 
		nodes_to_move = generate_combinations (list(G_primitive.nodes()), range(1, maxNodes+1))
		
		# choose blocks to move to 
		for nodescombo in nodes_to_move:
			# print('nodes to move', nodescombo)
			parts_to_move = [p for p in itertools.product(list(partDict.keys()), repeat=len(nodescombo))]   # permutation with repetition
			# parts_to_move = itertools.permutations(partDict.keys(), len(nodescombo))                      # permutation without repetition

			for partcombo in parts_to_move:
				valid = True 
				# print('parts to move to', partcombo)
				for idx, node_to_move in enumerate(nodescombo):
					node_part = get_part (partDict, node_to_move)
					# print(node_to_move, node_part, partcombo[idx])
					if node_part == partcombo[idx]:
						valid = False
						# print('invalid move')
				if valid: 
					# print('move is valid')
					partDict_tmp = copy.deepcopy(partDict)
					# print('original partdict', partDict_tmp)
					for idx, node_to_move in enumerate(nodescombo):
						node_part = get_part (partDict_tmp, node_to_move)
						partDict_tmp[node_part].remove(node_to_move)
						partDict_tmp[partcombo[idx]].append(node_to_move)
					# print('new partdict', partDict_tmp)
					part_sizes = [len(partDict_tmp[el]) for el in partDict_tmp]
					# print('part sizes', part_sizes)
					# check if all partitions are within size constraints after shifting
					if all(s <= Smax for s in part_sizes) and all(s >= Smin for s in part_sizes):
						T = max (calc_signal_path(G, innodes, outnodes, partDict_tmp))
						# print('highest T', T)
						
						C = cal_cut (G_primitive, partDict_tmp)
						
						# check if modules satisfy constraints
						part_opt_format = (C, [get_part(partDict_tmp, n) for n in G_primitive.nodes()])
						matrix, partG = partition_matrix (G, part_opt_format, innodes+outnodes, primitive_only)
						# print(matrix)
						loop_free = check_cycles(partG)
						# print(loop_free)
						motif_allowed = check_motif_allowed(matrix, motif_constraint)

						if loop_free and motif_allowed:
							# print('loop free and motif valid')
							solN += 1
							if T < minT:
								# print('optimized T from', minT, 'to', T)
								if os.path.exists(solfile):
									f_out = open(solfile, 'a')
								else:
									f_out = open(solfile, 'w')
								f_out.write('T\t'+str(T)+'\n')
								f_out.write('cut\t'+str(C)+'\n')
								for part in partDict_tmp.keys():
									f_out.write('Partition '+str(part)+'\t'+','.join(partDict_tmp[part])+'\n')
								f_out.close()
							else:
								if os.path.exists(potentialsolfile):
									f_out = open(potentialsolfile, 'a')
								else:
									f_out = open(potentialsolfile, 'w')
								f_out.write('T\t'+str(T)+'\n')
								f_out.write('cut\t'+str(C)+'\n')
								for part in partDict_tmp.keys():
									f_out.write('Partition '+str(part)+'\t'+','.join(partDict_tmp[part])+'\n')
								f_out.close()
	return solN

def optimize_signal_simulated_annealing (G, innodes, outnodes, primitive_only, cut, partDict, maxNodes, motif_constraint, outdir):
	""" 
	optimize based on signal travel time from inputs to the output 
	1. calculate the times that inputs have to traverse cell boundries 
	2. optmize the max traverse to be as low as possible 
	"""

	Smin, Smax = 3, 7
	nonprimitives = innodes + outnodes 

	G_primitive = copy.deepcopy(G)
	# if only considering primitives during partition and optimization
	if primitive_only: 
		for node in nonprimitives:
			if node in G.nodes():
				G_primitive.remove_node(node)
	print(G_primitive.nodes())
	# calculate original signal traverse time
	T = calc_signal_path (G, innodes, outnodes, partDict)
	minTC = max(T) + cut
	print(minTC)
	# get original partition matrix
	part_opt_format = (cut, [get_part(partDict, n) for n in G_primitive.nodes()])
	matrix, partG = partition_matrix (G, part_opt_format, nonprimitives, primitive_only)
	print('original partition', partDict)
	print('matrix of original partition', matrix)
	loop_free = check_cycles(partG)
	motif_allowed = check_motif_allowed(matrix, motif_constraint)
	# print('initial partition loop free', loop_free)
	# print('initial partition motif allowed', motif_allowed)

	## optimize traverse times 

	# make a directory to store results
	if os.path.exists(outdir):
		shutil.rmtree(outdir)
		os.mkdir(outdir)
	else: 
		os.mkdir(outdir)

	if len(G_primitive.nodes()) > Smax:

		# set parameters for simulated annealing
		Tmax = 2
		C    = 1e-4
		bestT_dict = dict()
		bestpartDict_all = dict()

		for i in range(1,2):  # for given number of trajectories
			print('iteration', i)
			bestpartDict = copy.deepcopy(partDict)
			# print(bestpartDict)
			bestT_list = [minTC]

			for t in range(1,10000):  # timestep in simulated annealing

				# at each timestep, choose a swap that satisfies the gate number constraints of each cell 
				size_constraint = False
				while size_constraint == False:

					# randomly choose a cell 
					cell = random.choice(list(partDict.keys()))
					# cell = 0
					# print('choosing cell', cell)
					# generate a subnetwork of this chosen cell and its neighboring cells
					neighbors = []
					for e in partG.edges():
						c1, c2 = e[0], e[1]
						if e[0] == cell: neighbors.append(e[1])
						elif e[1] == cell: neighbors.append(e[0])
					# print(list(set(neighbors)))
					subG_cells = list(set(neighbors)) + [cell]
					# print('subgraph cells', subG_cells)

					subG_nodes = [n for n in G_primitive.nodes() if get_part(bestpartDict, n) in subG_cells]
					# print('nodes in subgraph', subG_nodes)
					
					# choose 1 to n (maxNodes) nodes form this pair to swap 
					nodes_to_move = random.sample(subG_nodes, random.choice(np.arange(1,3)))
					# nodes_to_move = ['82']
					# print(nodes_to_move)

					partDict_tmp = copy.deepcopy(bestpartDict)

					# swap the selected nodes to other cells in this partition 
					for node in nodes_to_move:
						# print('move node', node)
						node_part = get_part(partDict_tmp, node)
						# print('original node part', node_part)
						new_part = node_part
						while new_part == node_part:
							new_part = random.choice(subG_cells)
						# print(new_part)
						partDict_tmp[node_part].remove(node)
						partDict_tmp[new_part].append(node)
					# print(partDict_tmp)

					# check if all cells are within size constrains after shifting
					part_sizes = [len(partDict_tmp[cell]) for cell in partDict_tmp]
					# print('part sizes', part_sizes)

					size_constraint = all(s <= Smax for s in part_sizes) and all(s >= Smin for s in part_sizes)
					# print('satisfy size constraints', size_constraint)

				# generate subgraph of the chosen cell and its neightbors
				# innodes_subG, outnodes_subG = classify_nodes (G_primitive, subG_nodes)
				# print(innodes_subG, outnodes_subG)
				# subG = nx.Graph.subgraph(G_primitive, subG_nodes)

				# calculate traverse time of the subgraph
				# T_sub_ori = max(calc_signal_path_subgraph(G_primitive, innodes_subG, outnodes_subG, subG_nodes, bestpartDict))
				# T_sub_new = max(calc_signal_path_subgraph(G_primitive, innodes_subG, outnodes_subG, subG_nodes, partDict_tmp))
				
				T_new = max(calc_signal_path(G, innodes, outnodes, partDict_tmp))

				# print('original T in subgraph', T_sub_ori)
				# print('new T', T_new)

				# check if modules in subgraph satisfy constraints
				cut_np = cal_cut (G, partDict_tmp)  # cut of new partition
				cut_bp = cal_cut (G, bestpartDict)  # cut of best partitionr result
				
				part_opt_format_new = (cut_np, [get_part(partDict_tmp, n) for n in G_primitive.nodes()])
				matrix_new, partG_new = partition_matrix(G, part_opt_format_new, nonprimitives, primitive_only)
				subG_matrix_new, subG_partG_new = get_part_matrix_subG (matrix_new, partG_new, subG_cells)
				# print('subG of new part', subG_matrix_new)
				part_opt_format_best = (cut_bp, [get_part(bestpartDict, n) for n in G_primitive.nodes()])
				matrix_best, partG_best = partition_matrix(G, part_opt_format_best, nonprimitives, primitive_only)
				subG_matrix_best, subG_partG_best = get_part_matrix_subG (matrix_best, partG_best, subG_cells)
				# print('subG of best part', subG_matrix_best)
				
				subG_new_loop_free = check_cycles(subG_partG_new)
				# print('subgraph loop free', subG_new_loop_free)
				subG_new_motif_allowed = check_motif_allowed(subG_matrix_new, motif_constraint)
				# print('subgraph motif allowed', subG_new_motif_allowed)

				subG_best_loop_free = check_cycles(subG_partG_best)
				# print('new subgraph loop free',subG_best_loop_free)
				subG_best_motif_allowed = check_motif_allowed(subG_matrix_best, motif_constraint)
				# print('new subgraph motif allowed', subG_best_motif_allowed)


				# part_opt_format_subG_ori = (cut_bestG, [get_part(bestpartDict, n) for n in subG_nodes])
				# # print(part_opt_format_subG_ori)
				# matrix_subG_ori, partG_subG_ori = partition_matrix(subG, part_opt_format_subG_ori, innodes_subG+outnodes_subG, primitive_only)
				# # print(matrix_subG_ori)
				# loop_free_subG_ori = check_cycles(partG_subG_ori)
				# # print(loop_free_subG_ori)
				# motif_allowed_subG_ori = check_motif_allowed(matrix_subG_ori, motif_constraint)

				# part_opt_format_subG = (cut_subG, [get_part(partDict_tmp, n) for n in subG_nodes])
				# # print(part_opt_format_subG)
				# matrix_subG, partG_subG = partition_matrix(subG, part_opt_format_subG, innodes_subG+outnodes_subG, primitive_only)
				# # print(matrix_subG)
				# loop_free_subG = check_cycles(partG_subG)
				# # print(loop_free_subG)
				# motif_allowed_subG = check_motif_allowed(matrix_subG, motif_constraint)


				# determine whether to accept or reject swap
				Ti = Tmax * (math.exp(-C*t))
				# print('Ti', Ti)
				deltaTC = (T_new+cut_np+1) - minTC

				if Ti > 0:
					try:
						P = math.exp(-(deltaTC/Ti))
					except OverflowError: 
						P = 2
					# print('temperature', Ti)
					# print('Probability', P)
				else:
					P = math.exp(-(deltaTC/0.01))
				# print('P', P)

				if subG_new_loop_free and subG_new_motif_allowed:
					print('new part loop free and motif valid')
					if subG_best_loop_free and subG_best_motif_allowed:
						print('original part loop free and motif valid')
						if T_new+cut_np < minTC:
							# accept swap 
							print('T improved, swap accepted')
							bestpartDict = partDict_tmp
							minTC = T_new+cut_np
						else: 
							# reject/accept with probability
							print('T NOT improved')
							if random.random() < P:
								bestpartDict = partDict_tmp
								minTC = T_new+cut_np
								print('swap accepted')
					else: 
						print('original part not loop free and motif valid')
						if T_new+cut_np <= minTC:
							# accept swap 
							print('T improved or equal')
							bestpartDict = partDict_tmp
							minTC = T_new+cut_np
						else: 
							# reject/accept with probability
							print('T not improved')
							if random.random() < P:
								bestpartDict = partDict_tmp
								minTC = T_new+cut_np
								print('swap accepted')
				else:
					print('new part NOT loop free and motif valid')
					if random.random() < P:
						bestpartDict = partDict_tmp
						minTC = T_new+cut_np
						print('swap accepted')


				bestT_list.append(minTC)

			# print('bestT', bestT_list)
			bestT_dict[i] = bestT_list 
			bestpartDict_all[i] = bestpartDict

		# write best partition result and best T 
		f_out = open(outdir + 'minTC.txt', 'w')
		f_out2 = open(outdir + 'optimized_solution_minTC.txt', 'w')
		f_out.write('iteration\tminT\n')
		for i in bestT_dict:
			f_out.write(str(i)+'\t'+','.join([str(T) for T in bestT_dict[i]])+'\n')
			cut = cal_cut (G_primitive, bestpartDict_all[i])
			T = max(calc_signal_path (G, innodes, outnodes, bestpartDict_all[i]))
			f_out2.write('path\t'+str(i)+'\n')
			f_out2.write('T\t'+str(T)+'\n')
			f_out2.write('cut\t'+str(cut)+'\n')
			for part in bestpartDict_all[i]:
				f_out2.write('Partition '+str(part)+'\t'+','.join(bestpartDict_all[i][part])+'\n')


def optimize_signal_subnetwork (G, innodes, outnodes, primitive_only, cut, partDict, maxNodes, motif_constraint, outdir):
	""" 
	optimize based on signal travel time from inputs to the output 
	1. calculate the times that inputs have to traverse cell boundries 
	2. optmize the max traverse to be as low as possible 
	"""

	Smin, Smax = 3, 7
	nonprimitives = innodes + outnodes 

	G_primitive = copy.deepcopy(G)
	# if only considering primitives during partition and optimization
	if primitive_only: 
		for node in nonprimitives:
			if node in G.nodes():
				G_primitive.remove_node(node)

	# calculate original signal traverse time
	T = calc_signal_path (G, innodes, outnodes, partDict)
	minTC = max(T) 

	# get original partition matrix
	part_opt_format = (cut, [get_part(partDict, n) for n in G_primitive.nodes()])
	matrix, partG = partition_matrix (G, part_opt_format, nonprimitives, primitive_only)
	# print('original partition', partDict)
	# print('matrix of original partition', matrix)
	loop_free = check_cycles(partG)
	motif_allowed = check_motif_allowed(matrix, motif_constraint)
	# print('initial partition loop free', loop_free)
	# print('initial partition motif allowed', motif_allowed)

	## optimize traverse times 

	# make a directory to store results
	if os.path.exists(outdir):
		shutil.rmtree(outdir)
		os.mkdir(outdir)
	else: 
		os.mkdir(outdir)

	if len(G_primitive.nodes()) > Smax:

		bestT_dict = dict()
		bestpartDict_all = dict()

		for i in range(1,11):  # for given number of trajectories
			print('iteration', i)
			bestpartDict = copy.deepcopy(partDict)
			bestpartG = copy.deepcopy(partG)
			# print(bestpartDict)
			bestT_list = [minTC]
			locked_nodes = []

			for t in range(1,10000):  # timestep in simulated annealing

				# at each timestep, choose a swap that satisfies the gate number constraints of each cell 
				size_constraint = False
				while size_constraint == False:
					# randomly choose a cell 
					cell = random.choice(list(bestpartDict.keys()))
					# cell = 0
					# print('choosing cell', cell)
					# generate a subnetwork of this chosen cell and its neighboring cells
					neighbors = []
					for e in bestpartG.edges():
						c1, c2 = e[0], e[1]
						if e[0] == cell: neighbors.append(e[1])
						elif e[1] == cell: neighbors.append(e[0])
					# print(list(set(neighbors)))
					subG_cells = list(set(neighbors)) + [cell]
					# print('subgraph cells', subG_cells)

					subG_nodes = list( set([n for n in G_primitive.nodes() if get_part(bestpartDict, n) in subG_cells]) - set(locked_nodes) )
					# print('nodes in subgraph', subG_nodes)
					
					# choose 1 to n (maxNodes) nodes form this pair to swap 
					nodes_to_move = random.sample(subG_nodes, random.choice(np.arange(1,3)))
					# nodes_to_move = ['82']
					# print(nodes_to_move)

					partDict_tmp = copy.deepcopy(bestpartDict)

					# swap the selected nodes to other cells in this partition 
					for node in nodes_to_move:
						# print('move node', node)
						node_part = get_part(partDict_tmp, node)
						# print('original node part', node_part)
						new_part = node_part
						while new_part == node_part:
							new_part = random.choice(subG_cells)
						# print(new_part)
						partDict_tmp[node_part].remove(node)
						partDict_tmp[new_part].append(node)
					# print(partDict_tmp)

					# check if all cells are within size constrains after shifting
					part_sizes = [len(partDict_tmp[cell]) for cell in partDict_tmp]
					# print('part sizes', part_sizes)

					size_constraint = all(s <= Smax for s in part_sizes) and all(s >= Smin for s in part_sizes)
					# print('satisfy size constraints', size_constraint)

				# generate subgraph of the chosen cell and its neightbors
				# innodes_subG, outnodes_subG = classify_nodes (G_primitive, subG_nodes)
				# print(innodes_subG, outnodes_subG)
				# subG = nx.Graph.subgraph(G_primitive, subG_nodes)

				# calculate traverse time of the subgraph
				# T_sub_ori = max(calc_signal_path_subgraph(G_primitive, innodes_subG, outnodes_subG, subG_nodes, bestpartDict))
				# T_sub_new = max(calc_signal_path_subgraph(G_primitive, innodes_subG, outnodes_subG, subG_nodes, partDict_tmp))
				
				T_new = max(calc_signal_path(G, innodes, outnodes, partDict_tmp))

				# print('original T in subgraph', T_sub_ori)
				# print('new T', T_new)

				# check if modules in subgraph satisfy constraints
				cut_np = cal_cut (G, partDict_tmp)  # cut of new partition
				cut_bp = cal_cut (G, bestpartDict)  # cut of best partitionr result
				
				part_opt_format_new = (cut_np, [get_part(partDict_tmp, n) for n in G_primitive.nodes()])
				matrix_new, partG_new = partition_matrix(G, part_opt_format_new, nonprimitives, primitive_only)
				subG_matrix_new, subG_partG_new = get_part_matrix_subG (matrix_new, partG_new, subG_cells)
				# print('subG of new part', subG_matrix_new)
				part_opt_format_best = (cut_bp, [get_part(bestpartDict, n) for n in G_primitive.nodes()])
				matrix_best, partG_best = partition_matrix(G, part_opt_format_best, nonprimitives, primitive_only)
				subG_matrix_best, subG_partG_best = get_part_matrix_subG (matrix_best, partG_best, subG_cells)
				# print('subG of best part', subG_matrix_best)
				
				subG_new_loop_free = check_cycles(subG_partG_new)
				# print('subgraph loop free', subG_new_loop_free)
				subG_new_motif_allowed = check_motif_allowed(subG_matrix_new, motif_constraint)
				# print('subgraph motif allowed', subG_new_motif_allowed)

				subG_best_loop_free = check_cycles(subG_partG_best)
				# print('new subgraph loop free',subG_best_loop_free)
				subG_best_motif_allowed = check_motif_allowed(subG_matrix_best, motif_constraint)

				if subG_new_loop_free and subG_new_motif_allowed:
					if subG_best_loop_free and subG_best_motif_allowed:
						if T_new < minTC:
							# accept swap 
							print('both part loop free and motif valid')
							print('T improved, swap accepted')
							bestpartDict = partDict_tmp
							minTC = T_new
							bestpartG = partG_new
							locked_nodes.extend (nodes_to_move)

					else: 
						print('original part not loop free and motif valid')
						# if T_new+cut_np <= minTC:
						# accept swap 
						print('T improved or equal')
						bestpartDict = partDict_tmp
						minTC = T_new
						bestpartG = partG_new
						locked_nodes.extend (nodes_to_move)

				bestT_list.append(minTC)

			# print('bestT', bestT_list)
			bestT_dict[i] = bestT_list 
			bestpartDict_all[i] = bestpartDict

		# write best partition result and best T 
		f_out = open(outdir + 'minT.txt', 'w')
		f_out2 = open(outdir + 'solution.txt', 'w')
		f_out.write('iteration\tminT\n')
		for i in bestT_dict:
			f_out.write(str(i)+'\t'+','.join([str(T) for T in bestT_dict[i]])+'\n')
			cut = cal_cut (G_primitive, bestpartDict_all[i])
			T = max(calc_signal_path (G, innodes, outnodes, bestpartDict_all[i]))
			f_out2.write('path\t'+str(i)+'\n')
			f_out2.write('T\t'+str(T)+'\n')
			f_out2.write('cut\t'+str(cut)+'\n')
			for part in bestpartDict_all[i]:
				f_out2.write('Partition '+str(part)+'\t'+','.join(bestpartDict_all[i][part])+'\n')



def plot_simulated_annealing_result (G, nonprimitives, primitive_only, outdir):
	""" plot partition result using simulated annealing """

	minT_file = outdir + 'minT.txt'
	opt_file = outdir + 'solution.txt'

	solDict = load_optimized_sol (opt_file)
	for sol in solDict:
		print('cut', solDict[sol]['cut'])
		print('T', solDict[sol]['T'])
		print(solDict[sol]['part'])

	G_primitive = copy.deepcopy(G)
	# if only considering primitives during partition and optimization
	if primitive_only: 
		for node in nonprimitives:
			if node in G.nodes():
				G_primitive.remove_node(node)

	for sol in [0]:
		part_opt_format = (solDict[sol]['cut'], [get_part(solDict[sol]['part'], n) for n in G_primitive.nodes()])
		print(solDict[sol]['part'])
		visualize_assignment_graphviz (G, outdir, part_opt_format, nonprimitives, primitive_only, outdir + 'part_opt_trajectory'+str(sol+1)+'.txt')

		# minT = min([solDict[s]['T'] for s in solDict])
		# minC = min([solDict[s]['cut'] for s in solDict if solDict[s]['T'] == minT])  # from all sols with minT, choose one with the lowest C
		# # print('minT', minT, 'minC', minC)
		# for s in solDict:
		# 	# print('solution', s)
		# 	# print(solDict[s])
		# 	if solDict[s]['T'] == minT and solDict[s]['cut'] == minC:
		# 		part_opt_format = (cut, [get_part(solDict[s]['part'], n) for n in G_primitive.nodes()])
		# 		record_opt_part (G, path, part_opt_format, solDict[s], nonprimitives, solN)
		# 		solN += 1

def optimize_signal_dfs (G, innodes, outnodes, primitive_only, cut, partDict, maxNodes, motif_constraint, outdir):

	Smin, Smax = 3, 7
	nonprimitives = innodes + outnodes 

	G_primitive = copy.deepcopy(G)
	# if only considering primitives during partition and optimization
	if primitive_only: 
		for node in nonprimitives:
			if node in G.nodes():
				G_primitive.remove_node(node)

	# calculate original signal traverse time
	T = calc_signal_path (G, innodes, outnodes, partDict)
	minT = max(T)

	# get original partition matrix
	part_opt_format = (cut, [get_part(partDict, n) for n in G_primitive.nodes()])
	matrix, partG = partition_matrix (G, part_opt_format, nonprimitives, primitive_only)
	print('original partition', partDict)
	print('matrix of original partition', matrix)
	loop_free = check_cycles(partG)
	motif_allowed = check_motif_allowed(matrix, motif_constraint)
	# print('initial partition loop free', loop_free)
	# print('initial partition motif allowed', motif_allowed)
	
	out_cell = []
	for cell in partG.nodes():
		if partG.out_degree(cell) == 0:
			out_cell.append(cell)
	print(out_cell)
	lp = 0
	for cell in partG.nodes():
		order = nx.dfs_tree(partG, source=cell)
		if len(order) > lp:
			lp = len(order)
		print(order.nodes())




def run_metis (G, nparts):
	# this function partitions a network into k subgraphs 
	partition = metis.part_graph(G, nparts = nparts, recursive=True)
	return partition

def find_best_n (G, nonprimitives, primitive_only): 
	""" for a given G, obtain the optimal nparts that gives the minimum cuts """
	# remove input and output nodes 
	G_primitive = copy.deepcopy(G)
	if primitive_only:
		for node in nonprimitives:
			if node in G_primitive.nodes():
				G_primitive.remove_node(node)

	if len(list(G_primitive.nodes())) > 7:
		min_nodes = 3
		max_nodes = 7
		nparts = list(range(int(len(G_primitive.nodes())/max_nodes)+1, int(len(G_primitive.nodes())/min_nodes)+1))
		# print('nparts', nparts)
		# find n (block number) that gives the minimum cuts 
		mincuts, part_opt = 10000, []
		for n in nparts:
			partition = run_metis (G_primitive, n)
			if partition[0] < mincuts:
				part_opt = partition 
				mincuts = part_opt[0]
	else: 
		part_opt = (0, [0]*len(list(G_primitive.nodes())))
	# print(part_opt)

	return part_opt

def partition_to_nparts (G, nonprimitives, primitive_only, npart):
	""" partition the circuit to n parts """
	# if primitive only, remove input and output nodes 
	G_primitive = copy.deepcopy(G)
	if primitive_only:
		for node in nonprimitives:
			if node in G_primitive.nodes():
				G_primitive.remove_node(node)

	part = run_metis (G_primitive, npart)
	return part


def partition_mincut_coi (dag_file, nonprimitives, primitive_only, outdir):
	""" partition dag such that cut is minimized """
	G = nx.Graph()
	G = nx.read_edgelist(dag_file, nodetype=str)

	part_opt = find_best_n(G, nonprimitives, primitive_only)

	# # write partition solution
	G_primitive = copy.deepcopy(G)
	for node in nonprimitives:
		if node in G_primitive.nodes():
			G_primitive.remove_node(node)

	outfile = outdir+'/part_solns.txt'

	f_out = open(outfile, 'w')
	f_out.write('cut\t'+str(part_opt[0])+'\n')

	# try:
	for part in range(max(part_opt[1])+1):
		nodeIdx = [a for a, b in enumerate(part_opt[1]) if b == part]
		nodes = [list(G.nodes())[n] for n in nodeIdx] 
		f_out.write('Partition '+str(part)+'\t')
		f_out.write(','.join([str(n) for n in nodes])+'\n')

	# visualize partition results
	outfig = outdir+'/part'
	outfile2 = outdir + '/part_matrix.npy'
	DiG = nx.read_edgelist(dag_file, nodetype=str, create_using=nx.DiGraph)
	visualize_assignment_graphviz (DiG, outdir, part_opt, nonprimitives, primitive_only, outfig)


def partition_nparts_coi (dag_file, nonprimitives, primitive_only, outdir):
	""" partition circuit into all possible nparts """
	G = nx.Graph()
	G = nx.read_edgelist(dag_file, nodetype=str)


	G_primitive = copy.deepcopy(G)
	if primitive_only:
		for node in nonprimitives:
			if node in G_primitive.nodes():
				G_primitive.remove_node(node)

	if len(list(G_primitive.nodes())) > 7:
		min_nodes = 3
		max_nodes = 7

		if os.path.exists(outdir+'/nparts'):
			shutil.rmtree(outdir+'/nparts')
			os.mkdir(outdir+'/nparts')
		else:
			os.mkdir(outdir+'/nparts')

		nparts = list(range(int(len(G_primitive.nodes())/max_nodes)+1, int(len(G_primitive.nodes())/min_nodes)+1))
		for n in nparts:
			part_opt = partition_to_nparts (G, nonprimitives, primitive_only, n)

			os.mkdir(outdir+'/nparts/'+str(n))
			outfile = outdir+'/nparts/'+str(n)+'/part_solns.txt'

			f_out = open(outfile, 'w')
			f_out.write('cut\t'+str(part_opt[0])+'\n')
                                                                                                                                      
			# try:
			for part in range(max(part_opt[1])+1):
				nodeIdx = [a for a, b in enumerate(part_opt[1]) if b == part]
				nodes = [list(G_primitive.nodes())[n] for n in nodeIdx] 
				f_out.write('Partition '+str(part)+'\t')
				f_out.write(','.join([str(n) for n in nodes])+'\n')

			# visualize partition results
			outfile = outdir+'/nparts/'+str(n)+'/part'
			outfile2 = outdir+'/nparts/'+str(n)+'/part_matrix.npy'
			DiG = nx.read_edgelist(dag_file, nodetype=str, create_using=nx.DiGraph)
			visualize_assignment_graphviz (DiG, outdir, part_opt, nonprimitives, primitive_only, outfile)


def visualize_graph_graphviz (G, nonprimitives, primitive_only, outfile):
	""" visualize DAG to be partitioned """
	fig = plt.figure(figsize=(8, 8))
	ax = fig.add_subplot(1,1,1)

	pos = graphviz_layout(G, prog='dot')
	nx.draw(G, pos, nodelist=list(G.nodes()), with_labels=True, node_color='blue')
	nx.draw(G, pos, nodelist=nonprimitives, node_color='lightgrey')
	plt.savefig(outfile+'/DAG.pdf', dpi=200)
	plt.show()



def visualize_assignment_graphviz (G, outdir, partition, nonprimitives, primitive_only, outfile): 
	""" visualize the partitioned graph with each color representing a partitioned block and shade representing toxicity"""
	fig = plt.figure(figsize=(16,8))

	# plot the partition assignment 
	ax = fig.add_subplot(1,2,1)
	# color = ['#bac964', '#438a5e', '#ffcbcb', '#f7d1ba', '#dbe3e5']
	color = sns.color_palette('hls', n_colors=max(partition[1])+1)

	G_primitive = copy.deepcopy(G)
	if primitive_only: 
		for node in nonprimitives:
			if node in G_primitive.nodes():
				G_primitive.remove_node(node)

	nx.nx_agraph.write_dot(G, outfile+'.dot')
	pos = graphviz_layout(G, prog='dot')
	for i in range(max(partition[1])+1):
		print('partition', i)
		nodeIdx = [a for a, b in enumerate(partition[1]) if b == i]
		nodes = [list(G_primitive.nodes())[n] for n in nodeIdx]
		print(nodes)
		# nx.draw(G, pos, nodelist=nodes, node_size=5, font_size=0.5, arrowsize=1, width=0.5, with_labels=True, node_color=color[i])
		nx.draw(G, pos, nodelist=nodes, with_labels=True, node_color=color[i])

	# nx.draw(G, pos, nodelist=nonprimitives, node_color='lightgrey')

	# plot partitioned cells
	ax2 = fig.add_subplot(1,2,2)

	matrix, partG =  partition_matrix (G, partition, nonprimitives, primitive_only)
	print(matrix)

	loops, nloops = [], []
	for e in partG.edges():
		if (e[1], e[0]) in partG.edges():
			loops.append(e)
		else: 
			nloops.append(e)
	pos = graphviz_layout(partG, prog='dot')
	nx.draw_networkx_nodes(partG, pos, node_size=300, node_color='#b18ea6')
	labels = {n:n for n in partG.nodes()}
	nx.draw_networkx_labels(partG, pos, labels)
	# draw loops 
	nx.draw_networkx_edges(partG, pos, loops, connectionstyle='arc3, rad=0.1')
	# draw non-loops 
	nx.draw_networkx_edges(partG, pos, nloops)
	ax2.axis('off')

	plt.savefig(outfile+'.pdf', dpi=200)
	# plt.show()

def partition_matrix (DiG, partition, nonprimitives, primitive_only):

	# generate adjency matrix
	numParts = max(partition[1])+1
	matrix = np.zeros(shape=(numParts, numParts))

	G_primitive = copy.deepcopy(DiG)
	if primitive_only:
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

def get_part_matrix_subG (matrix, partG, subG_cells):
	""" subset the matrix to only include cells in subgraph, and remove cells not from subgraph from partG """ 
	submatrix = np.take(matrix, subG_cells, axis=0)      # take rows of subgraph cells
	submatrix = np.take(submatrix, subG_cells, axis=1)   # take cols of subgraph cells

	partG_subG = copy.deepcopy(partG)
	for cell in partG:
		if cell not in subG_cells:
			partG_subG.remove_node(cell)

	return submatrix, partG_subG


def record_opt_part (G, path, part_opt, sol, nonprimitives, primitive_only, solN):

	# plot optimal partition graph

	visualize_assignment_graphviz (G, path, part_opt, nonprimitives, primitive_only, path+'/part_opt'+str(solN))

	# record optimal partition results
	f_out = open(path + '/part_opt'+str(solN)+'.txt', 'w')
	f_out.write('T\t'+str(sol['T'])+'\n')
	f_out.write('cut\t'+str(sol['cut'])+'\n')
	for part in range(max(sol['part'])+1):
		f_out.write('Partition '+str(part)+'\t')
		f_out.write(','.join([str(n) for n in sol['part'][part]])+'\n')

	# record partG matrix
	matrix, partG = partition_matrix (G, part_opt, nonprimitives, True)
	np.save(path + '/part_opt'+str(solN) + '.npy', matrix)

def record_ori_part (G, path, part_opt, partDict, nonprimitives, T, cut):
	visualize_assignment_graphviz (G, path, part_opt, nonprimitives, primitive_only, 'part_opt0')
	f_out = open(path+'/part_opt0.txt', 'w')
	f_out.write('T\t'+str(T)+'\n')
	f_out.write('cut\t'+str(cut)+'\n')
	for part in range(max(partDict)+1):
		f_out.write('Partition '+str(part)+'\t')
		f_out.write(','.join([str(n) for n in partDict[part]])+'\n')
	matrix, partG = partition_matrix (G, part_opt, nonprimitives, True)
	np.save(path + '/part_opt0.npy', matrix)


def get_optimal_results (G, G_primitive, cut, innodes, outnodes, primitive_only, partDict, constraint, path):
	signal_optimized_file = path + '/sigaltraverse_optimized_solution.txt'
	potential_sol_file = path + '/potential_solution.txt'
	nonprimitives = innodes+outnodes
	# check if the original partition fulfills the constraints
	T  = max( calc_signal_path (G, innodes, outnodes, partDict) )
	part_ori_format = (cut, [get_part(partDict, n) for n in G_primitive.nodes()])
	matrix, partG = partition_matrix (G, part_ori_format, nonprimitives, True)
	loop_free = check_cycles(partG)
	motif_allowed = check_motif_allowed(matrix, constraint)

	# print('original T', T)
	# print('cut', cut)
	# print('loop_free', loop_free)
	# print('motif', motif_allowed)

	solN = 1
	if os.path.exists(signal_optimized_file):
		print('optimized solution exists')
		solDict = load_optimized_sol (signal_optimized_file)
		minT = min([solDict[s]['T'] for s in solDict])
		minC = min([solDict[s]['cut'] for s in solDict if solDict[s]['T'] == minT])  # from all sols with minT, choose one with the lowest C
		# print('minT', minT, 'minC', minC)
		for s in solDict:
			# print('solution', s)
			# print(solDict[s])
			if solDict[s]['T'] == minT and solDict[s]['cut'] == minC:
				part_opt_format = (cut, [get_part(solDict[s]['part'], n) for n in G_primitive.nodes()])
				record_opt_part (G, path, part_opt_format, solDict[s], nonprimitives, primitive_only, solN)
				solN += 1

	else:
		if os.path.exists(potential_sol_file):
			# print('potential solutions')
			potsolDict = load_optimized_sol (potential_sol_file)
			minT = min([potsolDict[s]['T'] for s in potsolDict])
			minC = min([potsolDict[s]['cut'] for s in potsolDict if potsolDict[s]['T'] == minT])  # from all sols with minT, choose one with the lowest C
			
			if minT == T and minC < cut:	
				for s in potsolDict:
					# print('solution', s)
					# print(potsolDict[s])
					part_opt_format = (cut, [get_part(potsolDict[s]['part'], n) for n in G_primitive.nodes()])

					if potsolDict[s]['T'] == minT and potsolDict[s]['cut'] == minC:
						print('record potential solution that minimizes cut')
						record_opt_part (G, path, part_opt_format, potsolDict[s], nonprimitives, primitive_only, solN)
						solN += 1
			
			if solN == 1:
				print('no better solution than original')
				if motif_allowed and loop_free:
					print('recording original solution')
					record_ori_part (G, path, part_ori_format, partDict, nonprimitives, T, cut)
				else:  # if original sol doesn't satisfy constraints, pick one sol with lowest T and cut
					print('original solution does not satisfy constraints')
					for s in potsolDict:
						if potsolDict[s]['T'] == minT and potsolDict[s]['cut'] == minC:
							print('recording solution from potential solutions')
							part_opt_format = (cut, [get_part(potsolDict[s]['part'], n) for n in G_primitive.nodes()])
							# print(s, 'minT', minT, 'minC', minC)
							record_opt_part (G, path, part_opt_format, potsolDict[s], nonprimitives, primitive_only, solN)
							solN += 1


		else:
			# if original solution satisfies constraints
			if motif_allowed and loop_free:
				print('no potential solution file, recording original solution')
				record_ori_part (G, path, part_ori_format, partDict, nonprimitives, primitive_only, T, cut)
		

def check_all_solutions (graph, dag, nonprimitives, primitive_only, filepath):
	""" check all hc and lc solutions of npart graphs, to see how many are valid """
	f_out = open (filepath+ 'solution summary.txt', 'w')
	f_out.write('Graph	Nodes	Npart	Constraint	Iteration	Valid Motif_METIS	Valid Motif_Optimized	Cycle Free_METIS	Cycle Free_Optimized	T_Metis	T_Optimized	T_Delta\n')

	G_primitive = copy.deepcopy(dag)
	if primitive_only:
		for node in nonprimitives:
			if node in dag.nodes():
				G_primitive.remove_node(node)

	nparts = os.listdir(part_path)[1:]
	for npart in nparts:
		if npart.isdigit():
			print('npart', npart)
			nodes = len(list(dag.nodes()))
			for constraint in ['lc', 'hc']:
				if constraint == 'lc': motif_constraint = False
				else: motif_constraint = True
				
				# check if original partition satisfy constraints
				cut_o, part_o = load_part_sol (filepath+npart+'/part_solns.txt')
				part_opt_format_o = (cut_o, [get_part(part_o, n) for n in G_primitive.nodes()])
				matrix_o, partG_o = partition_matrix(dag, part_opt_format_o, nonprimitives, primitive_only)
				motif_allowed_o = check_motif_allowed(matrix_o, motif_constraint)
				loop_free_o = check_cycles(partG_o)

				TDict, cutDict, partDict = load_opt_part_sol (filepath+npart+'/optimized_'+constraint+'/solution.txt')
				minTDict = load_minT_sol (filepath+npart+'/optimized_'+constraint+'/minT.txt') 

				# check if motif_constraint is satisfied
				for iteration in TDict.keys():
					print('iteration', iteration)

					T_o = int(minTDict[iteration].split(',')[0])
					T = int(TDict[iteration])
					deltaT = T_o - T
					cut = int(cutDict[iteration])
					part = partDict[iteration]
					print(part)
					part_opt_format = (cut, [get_part(part, n) for n in G_primitive.nodes()])
					print(part_opt_format)
					for n in G_primitive.nodes():
						print(n, get_part(part, n))
					matrix, partG = partition_matrix (dag, part_opt_format, nonprimitives, primitive_only)
					# print('original partition', partDict)
					# print('matrix of original partition', matrix)
					loop_free = check_cycles(partG)
					motif_allowed = check_motif_allowed(matrix, motif_constraint)
					# compile results into one txt file 
					f_out.write('\t'.join([str(graph), str(nodes), npart, constraint])+'\t')
					f_out.write('\t'.join([str(iteration), str(motif_allowed_o), str(motif_allowed), str(loop_free_o), str(loop_free), str(T_o), str(T), str(deltaT)])+'\n')






if __name__ == '__main__':


	# load optimization improvement 
	# opt_results_file = '/Users/jgzhang/Work/Densmore lab/Partition/boolean circuit/optimization improvement - 4 input.txt'
	# opt_results = load_file (opt_results_file)
	# nonprimitives = ['a', 'b', 'c', 'd', 'e', 'out']

	# for i in range(60004,65530):
	# 	print(i)
	# 	path = '/Users/jgzhang/Programs/Cello2/sample-input/DNACompiler/4-input/'+str(i)
	# 	dag_file = path+'/DAG.edgelist'
	# 	# part_matrix = path+'/part_matrix.npy'
	# 	part_sol = path+'/part_solns.txt'
	# 	if os.path.exists(dag_file):
	# 		dag, innodes, outnodes = load_graph (dag_file)
	# 		nonprimitives = innodes + outnodes
	# 		# print('number of nodes', len(list(dag.nodes())))
	# 		# partG = load_part_graph (part_matrix)
	# 		G_primitive = copy.deepcopy(dag)
	# 		for node in nonprimitives:
	# 			if node in dag.nodes():
	# 				G_primitive.remove_node(node)

	# 		if os.stat(part_sol).st_size !=0 and len(list(dag.nodes())) <= 30:
	# 			cut, partDict = load_part_sol (part_sol)
			# 	part_opt_format = (cut, [get_part(partDict, n) for n in G_primitive.nodes()])
			# 	# # visualize_assignment_graphviz (dag, path, part_opt_format, nonprimitives, 'DAG')
			# 	T = calc_signal_path (dag, partDict)
			# 	# write_part_sol (path, partDict, T, cut)

			# 	if max(T) > 0:
			# 		solN = optimize_signal(dag, cut, partDict, primitive_only, 5, True, path+'/optimized_hc/')
			# 		if solN == 0:
			# 			solN_relaxed = optimize_signal(dag, cut, partDict, 5, False, path+'/optimized_lc/')

				# if os.path.exists(path+'/optimized_lc'):
				# 	print('optimized with low constraint')
				# 	get_optimal_results (dag, G_primitive, cut, innodes, outnodes, partDict, True, path+'/optimized_lc')
				# elif os.path.exists(path+'/optimized_hc'):
				# 	print('optimized with high constraint')
				# 	get_optimal_results (dag, G_primitive, cut, innodes, outnodes, partDict, False, path+'/optimized_hc')
				# else:	
				# 	get_optimal_results (dag, G_primitive, cut, innodes, outnodes, partDict, False, path)


				# if os.path.exists(path+'/optimized_lc'):
				# 	print('optimized with low constraint')

				# assign_qs (dag_file, part_matrix, part_sol)
				# optimize_cut_sa (dag, cut, partDict, [5, 0.5e-3], path)
				# numShifts = int( len(list(dag.nodes()))*0.2 )
				# print('number of shifts', numShifts)
				# optimize_cut_bf (dag, cut, numShifts, partDict)
			

	# circuit of interest
	# coi_path = '/Users/jgzhang/Programs/Cello2/sample-input/DNACompiler/bigcircuit/RCA4'
	coi_path = '/Users/jgzhang/Work/Densmore_lab/genetic-circuit-partitioning/2021.4/runs/benchmark/5-input-boolean-circuits'
	for graph in os.listdir(coi_path):
		if graph.startswith('.') == False:
			print(graph)
			outdir = coi_path + '/' + graph
			dag_file = outdir +'/DAG.edgelist'
			dag, innodes, outnodes = load_graph (dag_file) 
			nonprims = innodes+outnodes
			primitive_only = True
		# visualize_graph_graphviz (dag, nonprims, primitive_only, coi_path)
	# visualize_assignment_graphviz (dag, path, part_opt_format, nonprimitives, 'DAG')

	# G_primitive = copy.deepcopy(dag)
	# for node in nonprims:
	# 	if node in dag.nodes():
	# 		G_primitive.remove_node(node)

	# partition_coi (dag_file, nonprims, primitive_only, coi_path)
			# partition_nparts_coi (dag_file, nonprims, primitive_only, outdir)

			part_path = outdir+'/nparts/'
			for npart in os.listdir(part_path):
				if npart.isdigit():
					print('npart', npart)
					part_sol = part_path + npart+'/part_solns.txt'
					cut, partDict = load_part_sol (part_sol)
					# T = calc_signal_path (dag, innodes, outnodes, partDict)

					optimize_signal_subnetwork (dag, innodes, outnodes, True, cut, partDict, 5, False, part_path+npart+'/optimized_lc/')
					optimize_signal_subnetwork (dag, innodes, outnodes, True, cut, partDict, 5, True, part_path+npart+'/optimized_hc/')
				# 	# plot_simulated_annealing_result (dag, nonprims, True, part_path+npart+'/optimized_hc/')
				# 	# plot_simulated_annealing_result (dag, nonprims, True, part_path+npart+'/optimized_lc/')
				# check_all_solutions (graph, dag, nonprims, primitive_only, part_path)

	# # write_part_sol (part_path, partDict, T, cut)

	# # if max(T) > 0:
	# # 	solN = optimize_signal(dag, innodes, outnodes, cut, partDict, True, 5, True, part_path+'/optimized_hc/')
	# # 	if solN == 0:
	# # 		solN_relaxed = optimize_signal(dag, innodes, outnodes, cut, partDict, True, 5, False, part_path+'/optimized_lc/')

	# # if os.path.exists(part_path+'/optimized_lc'):
	# # 	print('optimized with low constraint')
	# # 	get_optimal_results (dag, G_primitive, cut, innodes, outnodes, True, partDict, True, part_path+'/optimized_lc')
	# # elif os.path.exists(part_path+'/optimized_hc'):
	# # 	print('optimized with high constraint')
	# # 	get_optimal_results (dag, G_primitive, cut, innodes, outnodes, True, partDict, False, part_path+'/optimized_hc')
	# # else:	
	# # 	get_optimal_results (dag, G_primitive, cut, innodes, outnodes, True, partDict, False, part_path)



	# dag_file = coi_path +'/DAG.edgelist'
	# dag, innodes, outnodes = load_graph (dag_file) 
	# nonprims = innodes+outnodes
	# primitive_only = True
	# part_path = coi_path+'/nparts/'
	# npart = '9'
	# part_sol = part_path + npart+'/part_solns.txt'
	# cut, partDict = load_part_sol (part_sol)
	# # T = calc_signal_path (dag, innodes, outnodes, partDict)

	# optimize_signal_subnetwork (dag, innodes, outnodes, True, cut, partDict, 5, False, part_path+npart+'/optimized_lc/')
	# optimize_signal_subnetwork (dag, innodes, outnodes, True, cut, partDict, 5, False, part_path+npart+'/optimized_hc/')
