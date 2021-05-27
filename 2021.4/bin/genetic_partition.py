#!/usr/bin/env python

#	Copyright (C) 2021 by
#	Jing Zhang <jgzhang@bu.edu>, Densmore Lab, Boston University
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.


# Load required modules
import csv
import random
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib
import networkx as nx
import metis
from collections import Counter
import numpy as np
import time
import numpy.linalg as la
import scipy.cluster.vq as vq
import itertools
import operator
import math
import copy
from mpmath import *
from itertools import chain
from itertools import product
from itertools import starmap
from functools import partial
import os
import seaborn as sns
import shutil
from networkx.drawing.nx_agraph import graphviz_layout
import re

##########################################
### create file names                  ###
##########################################

def edgelist_filename (settings, sample):
	return settings[sample]['graph_path']+'/DAG.edgelist'


##########################################
### load files                         ###
##########################################

def load_settings (filename):
	"""Load the settings file"""
	settings = {}
	data_reader = csv.reader(open(filename, 'rU'), delimiter='\t')
	# Ignore header
	header = next(data_reader)
	# Process each line
	for row in data_reader:
		if len(row) == len(header):
			sample = row[0]
			sample_data = {}
			for el_idx, el in enumerate(header[1:]):
				sample_data[el] = row[el_idx+1]
			settings[sample] = sample_data
	return settings


def load_graph (settings, sample):
	""" 
	read DAG edgelist, return DIRECTED graph, and input/output nodes 
	"""
	G = nx.read_edgelist (edgelist_filename (settings, sample), nodetype = str, create_using=nx.DiGraph())
	return G

def load_graph_undirected (settings, sample):
	"""
	read DAG edgelist, return UNDIRECTED graph, and input/output nodes 
	"""
	G = nx.Graph()
	G = nx.read_edgelist (edgelist_filename (settings, sample), nodetype=str)
	return G

def load_metis_part_sol (inputfile):
	"""
	read metis partition result 
	"""
	lines = [open(inputfile, 'r').read().strip("\n")][0].split('\n')
	cut = int( lines[0].split('\t')[1] )
	partDict = {}
	for line in lines[1:]: 
		tokens = line.split('\t')
		part = int( tokens[0].split(' ')[-1] )
		nodes = tokens[1].split(',')
		partDict[part] = nodes
	return cut, partDict

def load_opt_part_sol (inputfile):
	""" 
	load optimized solution 
	"""
	lines = [open(inputfile, 'r').read().strip("\n")][0].split('\n')
	iteration_idx = [lines.index(line) for line in lines if line.startswith('path') or line.startswith('sol')]
	solDict = {}
	partDict, TDict, cutDict = {}, {}, {}
	for i, idx in enumerate(iteration_idx):
		iteration = int (lines[idx].split('\t')[1] )
		solDict[iteration] = {}
		solDict[iteration]['T'] = int( lines[idx+1].split('\t')[1] )
		solDict[iteration]['cut'] = int( lines[idx+2].split('\t')[1] )
		solDict[iteration]['part'] = {}

		if idx!= iteration_idx[-1] and len(iteration_idx)!= 1: part_lines = lines[idx+3:iteration_idx[i+1]]
		else: part_lines = lines[idx+3:]

		for line in part_lines: 
			tokens = line.split('\t')
			part = int( tokens[0].split(' ')[-1] )
			nodes = tokens[1].split(',')
			solDict[iteration]['part'][part] = nodes
	return solDict


def load_minT_sol (inputfile):
	lines = [open(inputfile, 'r').read().strip("\n")][0].split('\n')
	minTDict = {}
	for line in lines[1:]:
		iteration, minT = int(line.split('\t')[0]), line.split('\t')[1]
		minTDict[iteration] = minT 
	return minTDict


###########################################
# Generate DAG edgelist from Cello's JSON output 
###########################################

def read_json(inputfile):
	lines = [open(inputfile, 'r').read().strip("\n")][0].split('\n')
	ports, gates = {}, {}
	for idx, line in enumerate(lines):
		line = line.strip()
		if line.startswith('"ports"'):
			p_s = idx
			searchlines = lines[idx+1:]
			for i, sl in enumerate(searchlines, idx):
				if sl.strip().startswith('"cells"'):
					p_e = i+1
		if line.startswith('"cells"'):
			g_s = idx
			searchlines = lines[idx+1:]
			for i, sl in enumerate(searchlines, idx):
				if sl.strip().startswith('"netnames"'):
					g_e = i
	# get information of inputs and outputs
	spacer = [idx+p_s+1 for idx, line in enumerate(lines[p_s+1:p_e]) if ': {' in line.strip()]
	for i, v in enumerate(spacer):
		# get names
		s = lines[v].strip()
		s = re.search('"(.*)"', s)
		el = s.group(1)
		ports[el] = {}
		# get directions 
		s = lines[v+1].strip()
		s = re.search('"direction": "(.*)"', s)
		direction = s.group(1)
		ports[el]['direction'] = direction
		# get bits
		s = lines[v+2].strip()
		bits = s.split('[')[1].split(']')[0].strip()
		ports[el]['bits'] = int(bits)
	# get information of gates 
	spacer = [idx+g_s+1 for idx, line in enumerate(lines[g_s+1:g_e]) if '$abc$' in line.strip()]
	for i, v in enumerate(spacer): 
		# get names
		s = int(lines[v].strip().split('"')[1].split('$')[-1])
		gates[s] = {}
		gates[s]['input'] = {}
		gates[s]['output'] = {}
		# search for attributes of this gate
		if i != len(spacer)-1:
			searchlines = lines[v:spacer[i+1]]
		else: 
			searchlines = lines[v:]
		for sl in searchlines:
			# get gate type
			if sl.strip().startswith('"type"'):
				gatetype = re.search('_(.*)_', sl.strip())
				gates[s]['type'] = gatetype.group(1)
			# get input(s)
			if sl.strip().startswith('"A": [') or sl.strip().startswith('"B": ['):
				port = re.search('"(.*)"', sl).group(1)
				bits = sl.split('[')[1].split(']')[0].strip()
				gates[s]['input'][port] = int(bits)
			# get output 
			if sl.strip().startswith('"Y": ['):
				port = re.search('"(.*)"', sl).group(1)
				bits = sl.split('[')[1].split(']')[0].strip()
				gates[s]['output'][port] = int(bits)
	return ports, gates

def synthesize_graph (ports, gates, outdir, t):
	G = nx.DiGraph()
	# start from the output, add edges
	edges = []

	for p in ports:
		if ports[p]['direction'] == 'output':
			b = ports[p]['bits']
			for g in gates:
				if b == gates[g]['output']['Y']:
					edges.append((g, p))
					# print('output', (g,p))

	for p in ports:
		if ports[p]['direction'] == 'input':
			b = ports[p]['bits']
			for g in gates:
				if b == gates[g]['input']['A']:
					edges.append((p, g))
					# print('input', (p,g))
				if gates[g]['type'] == 'NOR':
					if b == gates[g]['input']['B']:
						edges.append((p, g))
						# print('input', (p,g))

	for g in gates:
		op = gates[g]['output']['Y']
		for sg in gates: 
			if gates[sg]['type'] == 'NOT':
				gin = [gates[sg]['input']['A']]
			else:
				gin = [gates[sg]['input']['A'], gates[sg]['input']['B']]
			if op in gin:
				edges.append((g, sg))
				# print('internal', (g, sg))

	for e in edges:
		G.add_edge(*e)

	# write graph
	# nx.write_adjlist(G, outdir+'/DAG.adjlist')
	nx.write_edgelist(G, outdir+'/DAG.edgelist')


###########################################
# Supporting Functions
###########################################

def get_nonprimitive_nodes (G):
	"""
	Obtain nonprimitive nodes of a DAG
	input nodes (in_nodes) - in_degree is 0
	output nodes (out_nodes) - out_degree is 0
	"""
	in_nodes, out_nodes = [], []
	for node in G.nodes():
		indegree = G.in_degree(node)
		outdegree = G.out_degree(node)
		if outdegree == 0:
			out_nodes.append(node)
		if indegree == 0:
			in_nodes.append(node)
	nonprimitives = in_nodes + out_nodes 
	return in_nodes, out_nodes, nonprimitives 

def get_G_primitive (G, nonprimitives): 
	""" 
	if primitive only is True, remove input and output nodes 
	"""
	G_primitive = copy.deepcopy(G)
	for node in nonprimitives: 
		G_primitive.remove_node(node)
	return G_primitive

def get_part (partDict, node): 
	"""
	get the subgraph name that node is partitioned into 
	""" 
	for part, nodelist in partDict.items(): 
		if node in nodelist: 
			return part

def calc_signal_path (G, in_nodes, out_nodes, partDict):
	""" 
	count the number of boundaries that each input has to pass in order to reach the output 
	"""
	crosslist = []
	for inode in in_nodes:
		for outnode in out_nodes:
			for path in sorted(nx.all_simple_edge_paths(G, inode, outnode)):
				if (all(e in list(G.edges()) for e in path)): # valid path that makes a directed path from source to target
					cross = 0
					for e in path:
						if (any(n in in_nodes+out_nodes for n in e)) == False:
							if get_part(partDict, e[0]) != get_part(partDict, e[1]):
								cross += 1
					crosslist.append(cross)
	return crosslist

def calc_signal_path2 (partG):
	"""
	count the number of boundaries that each input has to pass in order to reach the output from partG graph
	"""
	crosslist = []
	in_nodes, out_nodes, nonprimitives = get_nonprimitive_nodes (partG)
	for in_node in in_nodes:
		for out_node in out_nodes:
			for path in sorted(nx.all_simple_edge_paths(partG, in_node, out_node)):
				if (all(e in list(partG.edges()) for e in path)): # valid path that makes a directed path from source to target
					cross = len(path) 
					crosslist.append(cross)
	return crosslist

def cal_cut (G, partDict):
	""" 
	calculate the cut in a new partition 
	"""
	cut = 0
	for edge in G.edges():
		part0 = get_part(partDict, edge[0])
		part1 = get_part(partDict, edge[1])
		if part0 != part1: 
			cut += 1
	return cut

def generate_combinations (n, rlist):
	""" from n choose r elements """
	combs = [list(itertools.combinations(n, r)) for r in rlist]
	combs = [item for sublist in combs for item in sublist]
	return combs

def check_cycles (partG):
	""" 
	check if partitioned cells contain cycles/loops 
	"""
	try:
		cycles = nx.find_cycle(partG)
		cycle_free = False
	except: 
		cycle_free = True
		pass 
	return cycle_free

def check_motif_allowed(matrix, motif_constraint):
	""" 
	check if all cells have allowed communications 
	""" 
	motif_allowed = True
	out_deg = np.sum(matrix, axis=1)
	in_deg = np.sum(matrix, axis=0)
	summed_deg = np.sum(matrix, axis=1)+np.sum(matrix, axis=0)

	if len(motif_constraint) == 2:    # high constraint 
		max_in, max_out = int(motif_constraint[0]), int(motif_constraint[1])
		if max(in_deg) > max_in: # sum by rows 
			motif_allowed = False
		if max(out_deg) > max_out: # sum by cols
			motif_allowed = False
		if max(summed_deg) > max_in + max_out:
			motif_allowed = False
	else:                             # low constraint 
		if max(summed_deg) > int(motif_constraint[0]):
			motif_allowed = False
	return motif_allowed



###########################################
# Perform initial graph partition
###########################################

def run_metis (G, n):
	""" 
	partitions a network into n subgraphs using metis
	"""
	part = metis.part_graph(G, nparts = n, recursive=True)
	return part

def partition_nparts_wrapper (G, n, outdir):
	""" 
	partition circuit into all possible nparts (wrapper function)
	"""
	# partition circuit into nparts 
	if len(list(G.nodes())) > 1:

		part_opt = run_metis (G, n)
		outfile = outdir +'/part_solns.txt'
		f_out = open(outfile, 'w')
		f_out.write('cut\t'+str(part_opt[0])+'\n')
                                                                                                                                  
		for part in range(max(part_opt[1])+1):
			nodeIdx = [a for a, b in enumerate(part_opt[1]) if b == part]
			nodes = [list(G.nodes())[node] for node in nodeIdx] 
			f_out.write('Partition '+str(part)+'\t')
			f_out.write(','.join([str(node) for node in nodes])+'\n')
	
	return part_opt

def partition_matrix (G, partition):
	"""
	generate subgraph adjency matrix, and DAG graph representing the subgraph network 
	"""
	# generate adjency matrix
	numParts = max(partition[1])+1
	matrix = np.zeros(shape=(numParts, numParts))

	for edge in G.edges():
		v1, v2 = edge[0], edge[1]
		part_v1 = partition[1][list(G.nodes()).index(v1)]
		part_v2 = partition[1][list(G.nodes()).index(v2)]
		if part_v1 != part_v2:
			matrix[part_v1][part_v2] += 1

	# generate DAG representing cell-cell communication from the adjency matrix
	rows, cols = np.where(matrix != 0)
	edges = zip(rows.tolist(), cols.tolist())
	partG = nx.DiGraph() 
	partG.add_edges_from(edges) 

	return matrix, partG

def get_part_matrix_subG (matrix, partG, subG_cells):
	""" 
	subset the matrix to only include cells in subgraph, and remove cells not from subgraph from partG 
	""" 
	submatrix = np.take(matrix, subG_cells, axis=0)      # take rows of subgraph cells
	submatrix = np.take(submatrix, subG_cells, axis=1)   # take cols of subgraph cells

	partG_subG = copy.deepcopy(partG)
	for cell in partG:
		if cell not in subG_cells:
			partG_subG.remove_node(cell)

	return submatrix, partG_subG

def rank_connectivity (G, primitive_only, outdir):
	"""
	ranks the initial npart networks by its connectivity
	"""
	if primitive_only:
		in_nodes, out_nodes, nonprimitives  = get_nonprimitive_nodes (G)
		G_primitive = get_G_primitive (G, nonprimitives)
	else:
		G_primitive   = copy.deepcopy (G)
		nonprimitives = []

	nparts = os.listdir(outdir)
	mean_part_degree = {}

	for npart in nparts: 
		if npart.isdigit():
			cut, partDict = load_metis_part_sol (outdir+npart+'/part_solns.txt')
			part_opt = (cut, [get_part(partDict, n) for n in G_primitive.nodes()])
			matrix, partG = partition_matrix(G_primitive, part_opt)
			sum_row = matrix.sum(axis=1)
			sum_col = matrix.sum(axis=0)
			mean_degree = np.median(sum_col + sum_row.T)
			mean_part_degree[npart] = mean_degree 

	return sorted(mean_part_degree.items(), key=lambda x: x[1])


###########################################
# Visualization 
###########################################

def visualize_assignment_graphviz (G, partition, nonprimitives, primitive_only, outdir, iteration): 
	""" 
	visualize the partitioned graph with each color representing a partitioned block
	"""
	if primitive_only:
		G_primitive = get_G_primitive (G, nonprimitives)
	else: 
		G_primitive = copy.deepcopy(G)

	fig = plt.figure(figsize=(16,8))

	# plot the partition assignment 
	ax = fig.add_subplot(1,2,1)
	# color = ['#bac964', '#438a5e', '#ffcbcb', '#f7d1ba', '#dbe3e5']
	color = sns.color_palette('hls', n_colors=max(partition[1])+1)
	# color=cm.rainbow(np.linspace(0,1,max(partition[1])+1))
	color = list(matplotlib.colors.cnames.values())

	nx.nx_agraph.write_dot(G, outdir+'/'+str(iteration)+'_DAG_part.dot')
	pos = graphviz_layout(G, prog='dot')
	for i in range(max(partition[1])+1):
		nodeIdx = [a for a, b in enumerate(partition[1]) if b == i]
		nodes = [list(G_primitive.nodes())[n] for n in nodeIdx]
		if len(list(G.nodes())) < 30: 
			nx.draw (G, pos, nodelist=nodes, with_labels=True, node_color=color[i])
		else:
			nx.draw (G, pos, nodelist=nodes, node_size=5, font_size=0.5, arrowsize=1, width=0.5, with_labels=True, node_color=color[i])

	# plot partitioned cells
	ax2 = fig.add_subplot(1,2,2)

	matrix, partG =  partition_matrix (G_primitive, partition)

	loops, nloops = [], []
	for e in partG.edges():
		if (e[1], e[0]) in partG.edges(): loops.append(e)
		else: nloops.append(e)

	pos = graphviz_layout(partG, prog='dot')
	nx.draw_networkx_nodes(partG, pos, node_size=300, node_color='#b18ea6')
	labels = {n:n for n in partG.nodes()}
	nx.draw_networkx_labels(partG, pos, labels)
	# draw loops 
	nx.draw_networkx_edges(partG, pos, loops, connectionstyle='arc3, rad=0.1')
	# draw non-loops 
	nx.draw_networkx_edges(partG, pos, nloops)
	ax2.axis('off')

	plt.savefig(outdir+'/'+str(iteration)+'_DAG_part.pdf', dpi=200)
	# plt.show()


###########################################
# Optimize partition result to satisfy motif constraints
###########################################


def optimize_signal_subnetwork (G, primitive_only, S_bounds, cut, partDict, maxNodes, motif_constraint, trajectories, outdir):
	""" 
	optimize based on signal travel time from inputs to the output 
	1. calculate the times that inputs have to traverse cell boundries 
	2. optmize the max traverse to be as low as possible 
	"""

	Smin, Smax = int (S_bounds[0]), int (S_bounds[1])

	if primitive_only:
		in_nodes, out_nodes, nonprimitives  = get_nonprimitive_nodes (G)
		G_primitive = get_G_primitive (G, nonprimitives)
	else:
		G_primitive   = copy.deepcopy (G)
		nonprimitives = []

	# calculate original signal traverse time
	T = calc_signal_path (G, in_nodes, out_nodes, partDict)
	minT = max(T) 
	# print(minT)
	# get original partition matrix
	part_opt = (cut, [get_part(partDict, n) for n in G_primitive.nodes()])
	matrix, partG = partition_matrix (G_primitive, part_opt)

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

		for i in range(1, trajectories+1):  # for given number of trajectories
			print('iteration', i)
			bestpartDict = copy.deepcopy(partDict)
			bestpartG    = copy.deepcopy(partG)
			bestMatrix   = copy.deepcopy(matrix)
			bestT_list = [minT]
			locked_nodes = []

			for t in range(1,10000):  # timestep
				# at each timestep, choose a swap that satisfies the gate number constraints of each cell 
				# print('t', t)
				# get partG of the best part
				cut_bp = cal_cut (G, bestpartDict)
				part_opt_format_best = (cut_bp, [get_part(bestpartDict, n) for n in G_primitive.nodes()])
				matrix_best, partG_best = partition_matrix(G_primitive, part_opt_format_best)
				
				# get subnetworks that do not satisfy constraints
				cell_unmet_const = []
				cell_met_const   = []
				for cell in list(bestpartG.nodes()):
					neighbors = []
					for e in list(bestpartG.edges()):
						c1, c2 = e[0], e[1]
						if e[0] == cell: neighbors.append(e[1])
						elif e[1] == cell: neighbors.append(e[0])
					subG_cells = list(set(neighbors)) + [cell]
					subG_matrix_best, subG_partG_best = get_part_matrix_subG (matrix_best, partG_best, subG_cells)
					subG_best_loop_free = check_cycles(subG_partG_best)
					subG_best_motif_allowed = check_motif_allowed(subG_matrix_best, motif_constraint)
					if subG_best_loop_free and subG_best_motif_allowed: 
						cell_met_const.append (cell)
					else: 
						cell_unmet_const.append (cell)
				# print('constraint ok', cell_met_const)
				# print('constraint not ok', cell_unmet_const)
				size_constraint = False
				while size_constraint == False:
					# randomly choose a cell 
					try:
						if random.uniform(0, 1) < 0.2: 
							cell = random.choice(cell_met_const)
						else: 
							cell = random.choice(cell_unmet_const)
					except IndexError: 
						cell = random.choice(list(bestpartG.nodes()))
					# print('choosing cell', cell)
					# generate a subnetwork of this chosen cell and its neighboring cells
					neighbors = []
					for e in list(bestpartG.edges()):
						c1, c2 = e[0], e[1]
						if e[0] == cell: neighbors.append(e[1])
						elif e[1] == cell: neighbors.append(e[0])
					# print(neighbors)
					subG_cells = list(set(neighbors)) + [cell]
					# print('subgraph cells', subG_cells)

					subG_nodes = list( set([n for n in G_primitive.nodes() if get_part(bestpartDict, n) in subG_cells]) - set(locked_nodes) )
					# print('nodes in subgraph', subG_nodes)
					
					# choose 1 to n (maxNodes) nodes form this pair to swap 
					trial, have_nodes_to_move = 0, False
					while  have_nodes_to_move == False:
						try: 
							nodes_to_move = random.sample(subG_nodes, random.choice(np.arange(1,maxNodes+1)))
							have_nodes_to_move = True
						except ValueError: 
							have_nodes_to_move = False
							trial += 1
							if trial > 50: break      

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

					size_constraint = all(s <= Smax for s in part_sizes) and all(s >= Smin for s in part_sizes)
					# print('size constraint', size_constraint)
				# T_new = max(calc_signal_path(G, in_nodes, out_nodes, partDict_tmp))

				# print('original T in subgraph', T_sub_ori)
				# print('new T', T_new)

				# check if modules in subgraph satisfy constraints
				cut_np = cal_cut (G, partDict_tmp)  # cut of new partition
				cut_bp = cal_cut (G, bestpartDict)  # cut of best partitionr result
				
				part_opt_format_new = (cut_np, [get_part(partDict_tmp, n) for n in G_primitive.nodes()])
				matrix_new, partG_new = partition_matrix(G_primitive, part_opt_format_new)
				try:
					T_new = max(calc_signal_path2 (partG_new))				
					subG_matrix_new, subG_partG_new = get_part_matrix_subG (matrix_new, partG_new, subG_cells)
					# print('subG of new part', subG_matrix_new)
					subG_new_loop_free = check_cycles(subG_partG_new)
					# print('subgraph loop free', subG_new_loop_free)
					subG_new_motif_allowed = check_motif_allowed(subG_matrix_new, motif_constraint)
					# print('subgraph motif allowed', subG_new_motif_allowed)
					subG_matrix_best, subG_partG_best = get_part_matrix_subG (matrix_best, partG_best, subG_cells)
					subG_best_loop_free = check_cycles(subG_partG_best)
					subG_best_motif_allowed = check_motif_allowed(subG_matrix_best, motif_constraint)

					if subG_new_loop_free and subG_new_motif_allowed:
						if subG_best_loop_free and subG_best_motif_allowed:
							if T_new < minT:
								# accept swap 
								print('both part loop free and motif valid')
								print('T improved, swap accepted')
								bestpartDict = partDict_tmp
								bestpartG    = partG_new
								bestMatrix   = matrix_new
								minT = T_new
								locked_nodes.extend (nodes_to_move)

						else: 
							print('original part not loop free and motif valid')
							# if T_new+cut_np <= minT:
							# accept swap 
							print('T improved or equal')
							bestpartDict = partDict_tmp
							bestpartG    = partG_new
							bestMatrix   = matrix_new
							minT = T_new
							locked_nodes.extend (nodes_to_move)

					bestT_list.append(minT)
				except ValueError: 
					pass
			# print('bestT', bestT_list)
			bestT_dict[i] = bestT_list 
			bestpartDict_all[i] = bestpartDict

		print('recording solution')
		# write best partition result and best T 
		f_out = open(outdir + 'minT.txt', 'w')
		f_out2 = open(outdir + 'part_solns.txt', 'w')
		f_out.write('iteration\tminT\n')
		for i in bestT_dict:
			f_out.write(str(i)+'\t'+','.join([str(T) for T in bestT_dict[i]])+'\n')
			cut = cal_cut (G_primitive, bestpartDict_all[i])
			T = max(calc_signal_path (G, in_nodes, out_nodes, bestpartDict_all[i]))
			f_out2.write('path\t'+str(i)+'\n')
			f_out2.write('T\t'+str(T)+'\n')
			f_out2.write('cut\t'+str(cut)+'\n')
			for part in bestpartDict_all[i]:
				f_out2.write('Partition '+str(part)+'\t'+','.join(bestpartDict_all[i][part])+'\n')

def check_all_solutions (graph, G, nonprimitives, primitive_only, high_constraint, low_constraint, outdir):
	""" 
	check all hc and lc solutions of npart graphs, to see how many are valid 
	"""
	f_out = open (outdir + 'solution summary.txt', 'w')
	f_out.write('\t'.join(['Graph', 'Nodes', 'Npart', 'Constraint', 'Iteration','Valid Motif_METIS', 'Valid Motif_Optimized', 'Cycle Free_METIS', 'Cycle Free_Optimized', 'T_Metis', 'T_Optimized', 'T_Delta'])+'\n')

	if primitive_only:
		in_nodes, out_nodes, nonprimitives  = get_nonprimitive_nodes (G)
		G_primitive = get_G_primitive (G, nonprimitives)
	else:
		G_primitive   = copy.deepcopy (G)
		nonprimitives = []

	# check if original partition satisfy constraints
	cut_o, part_o = load_metis_part_sol (outdir+npart+'/part_solns.txt')
	part_opt_format_o = (cut_o, [get_part(part_o, n) for n in G_primitive.nodes()])
	matrix_o, partG_o = partition_matrix(G_primitive, part_opt_format_o)
	motif_allowed_o = check_motif_allowed(matrix_o, motif_constraint)
	loop_free_o = check_cycles(partG_o)

	nparts = os.listdir(part_path)[1:]
	for npart in nparts:
		if npart.isdigit():
			print('npart', npart)
			nodes = len(list(G.nodes()))
			for constraint in ['lc', 'hc']:
				if constraint == 'lc': motif_constraint = low_constraint
				else: motif_constraint = high_constraint

				solDict = load_opt_part_sol (outdir+npart+'/optimized_'+constraint+'/solution.txt')
				minTDict = load_minT_sol (outdir+npart+'/optimized_'+constraint+'/minT.txt') 

				# check if motif_constraint is satisfied
				for iteration in solDict.keys():
					print('iteration', iteration)

					T_o = int(minTDict[iteration].split(',')[0])
					T = int(solDict[iteration]['T'])
					deltaT = T_o - T
					cut = int(solDict[iteration]['cut'])
					part = solDict[iteration]['part']
					print(part)
					part_opt_format = (cut, [get_part(part, n) for n in G_primitive.nodes()])
					print(part_opt_format)
					for n in G_primitive.nodes():
						print(n, get_part(part, n))
					matrix, partG = partition_matrix (G_primitive, part_opt_format_new)
					# print('original partition', partDict)
					# print('matrix of original partition', matrix)
					loop_free = check_cycles(partG)
					motif_allowed = check_motif_allowed(matrix, motif_constraint)
					# compile results into one txt file 
					f_out.write('\t'.join([str(graph), str(nodes), npart, constraint])+'\t')
					f_out.write('\t'.join([str(iteration), str(motif_allowed_o), str(motif_allowed), str(loop_free_o), str(loop_free), str(T_o), str(T), str(deltaT)])+'\n')

   
def optimize_signal_bruteforce (G, primitive_only, S_bounds, cut, partDict, maxNodes, motif_constraint, outdir):
	""" 
	optimize based on signal travel time from inputs to the output, search for all posible combinations
	1. calculate the times that inputs have to traverse cell boundries 
	2. optmize the max traverse to be as low as possible 
	"""

	Smin, Smax = int (S_bounds[0]), int (S_bounds[1])

	if primitive_only:
		in_nodes, out_nodes, nonprimitives  = get_nonprimitive_nodes (G)
		G_primitive = get_G_primitive (G, nonprimitives)
	else:
		G_primitive   = copy.deepcopy (G)
		nonprimitives = []

	# calculate original signal traverse time
	T = calc_signal_path (G, in_nodes, out_nodes, partDict)
	minT = max(T) 
	# get original partition matrix
	part_opt_o = (cut, [get_part(partDict, n) for n in G_primitive.nodes()])
	matrix_o, partG_o = partition_matrix (G_primitive, part_opt_o)
	loop_free_o = check_cycles(partG_o)
	motif_allowed_o = check_motif_allowed(matrix_o, motif_constraint)

	# make a directory to store results
	if os.path.exists(outdir):
		shutil.rmtree(outdir)
		os.mkdir(outdir)
	else: 
		os.mkdir(outdir)

	if minT > 0:
		solN = 0
		# store solutions 
		f_out = open(outdir + 'part_solns.txt', 'w')

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
					partDict_tmp = copy.deepcopy(partDict)
					# print('original partdict', partDict_tmp)
					for idx, node_to_move in enumerate(nodescombo):
						node_part = get_part (partDict_tmp, node_to_move)
						partDict_tmp[node_part].remove(node_to_move)
						partDict_tmp[partcombo[idx]].append(node_to_move)
					# print('new partdict', partDict_tmp)
					part_sizes = [len(partDict_tmp[el]) for el in partDict_tmp]
					# check if all partitions are within size constraints after shifting
					if all(s <= Smax for s in part_sizes) and all(s >= Smin for s in part_sizes):
						T = max (calc_signal_path(G, in_nodes, out_nodes, partDict_tmp))

						C = cal_cut (G_primitive, partDict_tmp)
						
						# check if modules satisfy constraints
						part_opt_format = (C, [get_part(partDict_tmp, n) for n in G_primitive.nodes()])
						matrix, partG = partition_matrix (G_primitive, part_opt_format)
						# print(matrix)
						loop_free = check_cycles(partG)
						# print(loop_free)
						motif_allowed = check_motif_allowed(matrix, motif_constraint)

						if loop_free and motif_allowed:
							if T <= minT: 
								solN += 1
								f_out.write('sol\t'+str(solN)+'\n'+'T\t'+str(T)+'\n'+'cut\t'+str(C)+'\n')
								for part in partDict_tmp.keys(): f_out.write('Partition '+str(part)+'\t'+','.join(partDict_tmp[part])+'\n')
							else: 
								if not (motif_allowed_o and loop_free_o):
									solN += 1
									f_out.write('sol\t'+str(solN)+'\n'+'T\t'+str(T)+'\n'+'cut\t'+str(C)+'\n')
									for part in partDict_tmp.keys(): f_out.write('Partition '+str(part)+'\t'+','.join(partDict_tmp[part])+'\n')
								

def determine_best_solution (G, primitive_only, high_constraint, low_constraint, outdir):
	""" 
	from all solutions, choose the one(s) that satisfies motif constraints (prioritize high constraints, then low), 
	while minimizing T (and cut)
	"""
	f_out = open (outdir + 'best_solns.txt', 'w')
	f_out.write('\t'.join(['Npart', 'Sol', 'Nodes', 'Constraint', 'Valid Motif_METIS', 'Valid Motif_Optimized', 'Cycle Free_METIS', 'Cycle Free_Optimized', 'T_Metis', 'T_Optimized', 'cut_Metis', 'cut_Optimized'])+'\n')

	if primitive_only:
		in_nodes, out_nodes, nonprimitives  = get_nonprimitive_nodes (G)
		G_primitive = get_G_primitive (G, nonprimitives)
	else:
		G_primitive   = copy.deepcopy (G)
		nonprimitives = []

	nparts = os.listdir(outdir)

	for constraint in ['lc', 'hc']:
		if constraint == 'lc': motif_constraint = low_constraint
		else: motif_constraint = high_constraint
		for npart in nparts:
			if npart.isdigit():
				if os.path.exists(outdir+npart+'/optimized_'+constraint+'/part_solns.txt'):
					# check if original partition satisfy constraints
					cut_o, part_o = load_metis_part_sol (outdir+npart+'/part_solns.txt')
					T_o = max(calc_signal_path (G, in_nodes, out_nodes, part_o))
					part_opt_o = (cut_o, [get_part(part_o, n) for n in G_primitive.nodes()])
					matrix_o, partG_o = partition_matrix(G_primitive, part_opt_o)
					motif_allowed_o = check_motif_allowed(matrix_o, motif_constraint)
					loop_free_o = check_cycles(partG_o)
					best_soln = [('0', part_opt_o)]
					minT      = T_o
					# print(best_soln)
					# load optimized solution
					solDict = load_opt_part_sol (outdir+npart+'/optimized_'+constraint+'/part_solns.txt')

					# check if motif_constraint is satisfied
					for iteration in solDict.keys():
						T = int(solDict[iteration]['T'])
						cut = int(solDict[iteration]['cut'])
						part = solDict[iteration]['part']
						part_opt = (cut, [get_part(part, n) for n in G_primitive.nodes()])

						matrix, partG = partition_matrix (G_primitive, part_opt)
						loop_free = check_cycles(partG)
						motif_allowed = check_motif_allowed(matrix, motif_constraint)

						# best soln so far
						matrix_bs, partG_bs = partition_matrix (G_primitive, best_soln[0][1])
						loop_free_bs = check_cycles(partG_bs)
						motif_allowed_bs = check_motif_allowed(matrix_bs, motif_constraint)

						if not (motif_allowed_bs and loop_free_bs):  # if best solution doesn't satisfy constraint
							if motif_allowed and loop_free:          # if new part does
								best_soln = [(iteration, part_opt)]
								# print(best_soln)
								minT      = T
						else:                                        # if best solution satisfies cnstraint
							if motif_allowed and loop_free:          # if new part does
								if T < minT:
									best_soln = [(iteration, part_opt)]
									minT      = T
								elif T == minT:
									if cut < best_soln[0][1][0]: 
										best_soln = [(iteration, part_opt)]
									elif cut == best_soln[0][1][0]:
										best_soln.append((iteration, part_opt)) 
										# print(best_soln)

					for soln in best_soln:
						# compile results
						matrix_bs, partG_bs = partition_matrix (G_primitive, soln[1])
						loop_free_bs = check_cycles(partG_bs)
						motif_allowed_bs = check_motif_allowed(matrix_bs, motif_constraint)
						if loop_free_bs and motif_allowed_bs:
							f_out.write('\t'.join([str(npart), str(soln[0]), str(len(list(G.nodes()))), constraint, str(motif_allowed_o), str(motif_allowed_bs), str(loop_free_o), str(loop_free_bs), str(T_o), str(minT), str(cut_o), str(soln[1][0])])+'\n')
							# visualize best solutions
							visualize_assignment_graphviz (G, soln[1], nonprimitives, primitive_only, outdir+npart+'/optimized_'+constraint, soln[0])



