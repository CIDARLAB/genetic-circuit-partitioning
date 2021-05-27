#!/usr/bin/env python

#	Copyright (C) 2020 by
#	Jing Zhang <jgzhang@bu.edu>, Densmore Lab, Boston University
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

# Required modules
import random
import matplotlib.pyplot as plt
from matplotlib import cm
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


##########################################
### create file names                  ###
##########################################

def edgelist_filename (settings, sample):
	return settings[sample]['graph_path']+'/DAG'+settings[sample]['n']+'.'+settings[sample]['k']+'.edgelist'

def DAG_edgelist_filename (settings, sample):
	return settings[sample]['graph_path']+'/DAG'+settings[sample]['n']+'.'+settings[sample]['k']+'.edgelist'

def gate_lib_filename (settings, sample):
	return settings[sample]['lib_path']+'gatelib.txt'

def sensor_lib_filename (settings, sample):
	return settings[sample]['lib_path']+'sensor.txt'

def logic_filename (settings, sample):
	return settings[sample]['lib_path']+'statelogic.txt'

def part_sol_filename (settings, sample):
	return settings[sample]['output_path']+'/DAG'+settings[sample]['n']+'.'+settings[sample]['k']+'/part_solutions.txt'

def part_outfig (settings, sample):
	return settings[sample]['output_path']+'/DAG'+settings[sample]['n']+'.'+settings[sample]['k']+'/part_fig.pdf'

def S_trajectory_filename (settings, sample):
	return settings[sample]['output_path']+'/DAG'+settings[sample]['n']+'.'+settings[sample]['k']+'/S_trajectory.txt'

def S_trajectory_plotname (settings, sample):
	return settings[sample]['output_path']+'/DAG'+settings[sample]['n']+'.'+settings[sample]['k']+'/S_trajectory.pdf'

def opt_gate_assignment_filename (settings, sample, run):
	return settings[sample]['output_path'+'/DAG'+settings[sample]['n']+'.'+settings[sample]['k']+'/optimal gate assignments_trajectory_'+str(run)+'.txt'



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

def load_gate_info (settings, sample):
	lines = [open(gate_lib_filename (settings, sample), 'r').read().strip("\n")][0].split('\n')
	gateLib = {}
	header = lines[0].split('\t')
	for line in lines[1:]: 
		tokens = line.split('\t')
		gate = tokens[0]
		gateLib[gate] = {}
		for idx, token in enumerate(tokens, 1):
			gateLib[gate][header[idx-1]] = token
	return gateLib

def load_sensor_info (settings, sample):
	lines = [open(sensor_lib_filename (settings, sample), 'r').read().strip("\n")][0].split('\n')
	sensorLib = {}
	for line in lines[1:]:
		tokens = line.split('\t')
		sensor = tokens[0]
		REU_off = float(tokens[1])
		REU_on = float(tokens[2])
		sensorLib[sensor] = {'on': REU_on, 'off': REU_off}
	return sensorLib

def load_logic (settings, sample):
	lines = [open(logic_filename (settings, sample), 'r').read().strip("\n")][0].split('\n')
	logic = {'PTac':[], 'PTet':[], 'PBAD':[]}
	for line in lines[1:]:0
		tokens = line.split('\t')[1].split(',')
		logic['PTac'].append(int(tokens[0]))
		logic['PTet'].append(int(tokens[1]))
		logic['PBAD'].append(int(tokens[2]))
	return logic

def read_edge_direction (settings, sample):
	""" this function reads the .edgelist file and returns the direction of G edges """
	lines = [open(DAG_edgelist_filename (settings, sample), 'r').read().strip("\n")][0].split('\n')
	edge_direct = []
	for line in lines: 
		words = line.split(' ')
		edge_direct.append([words[0], words[1]])
	return edge_direct

def load_gate_assignment (settings, sample, G):
	""" load partitioned circuit with biological gates assigned """
	lines = [open(S_trajectory_filename (settings, sample), 'r').read().strip("\n")][0].split('\n')
	# get the best run with the highest circuit score
	bestrun, bestS = 'NA', 0
	for line in lines:
		run, SList = line.split('\t')[0], [float(s) for s in line.split('\t')[1].split(',')]
		if SList[-1] > bestS:
			bestrun = run
	# load the gate assignment data from the best run
	lines = [open(opt_gate_assignment_filename (settings, sample, bestrun), 'r').read().strip("\n")][0].split('\n')

	# update G with gate assignment information 
	for line in lines[1:]:
		tokens = line.split('\t')
		node = tokens[0]
		G.nodes[node]['type'] = tokens[1]
		G.nodes[node]['gate'] = tokens[2]
		G.nodes[node]['params'] = tokens[3]
		G.nodes[node]['part'] = tokens[4]
	return G
 

##########################################
### define functions                   ###
##########################################

def sigmoid( ymin, ymax, Kd, n, x): 
	return ymin + (ymax - ymin)/(1 + math.pow( (x/Kd), n) )

def sigmoid_curve(xlist, ymin, ymax, Kd, n):
	ylist = [sigmoid(ymin, ymax, Kd, n, x) for x in xlist]
	return ylist



###########################################
# Generate random graphs                ###
###########################################

def generate_degree_seq(n):
	""" generate a degree sequence """
	# for a given number n, randomly generate n1 number of NOT (degree: 2), 
	# and n2 number of NOR (degree: 3) gates
	seq = [1]
	while sum(seq) %2 != 0:
		seq = [random.randint(2,3) for x in range(n)]
	return seq

def generate_degree_seq2(n):
	""" generate a degree sequence following prior knowledge of degree distribution """
	# as we know that the ratio of degree 3 nodes and degree 2 nodes is 2.16:1 in order to get a connected graph 
	# we can generate random degree seq satisfying this ratio to get successful graphs, also lowers the search space when n is large
	seq = [1]
	choice_list = [2]*32 + [3]*68
	while sum(seq) %2 != 0:
		seq = random.sample(choice_list, n)
	return seq

def generate_random_graph(n):
	""" for a given vertice number n, generate a connected graph """
	n0 = 4                                                 # define the number of primitives with degree 1 
	if n <= 50: 
		z = [1]*n0 + sorted(generate_degree_seq (n-n0))
	else: 
		z = [1]*n0 + sorted(generate_degree_seq2 (n-n0))    
	G = nx.configuration_model(z, create_using=nx.Graph)   # configuration model
	while nx.is_connected(G) == False:                     # make sure the graph is connected
		if n <= 50: 
			z = [1]*n0 + sorted(generate_degree_seq (n-n0))
		else: 
			z = [1]*n0 + sorted(generate_degree_seq2 (n-n0))
		G = nx.configuration_model(z, create_using=nx.Graph)
	return G, Counter(z)


def generate_random_graph (settings, sample): 
	""" 
	generate a connected graph, which follows these constraints -
	1. does not contain self loop or parallel edges 
	2. two degree 2 nodes are not connected to each other
	3. the graph is connected
	"""
    n_nonprimitives = len( settings[sample]['in_nodes'] ) + len( settings[sample]['out_nodes'] )
    n_primitives = int( settings[sample]['n'] ) - n_nonprimitives

	# generate a random graph in which two d2 nodes are not connected to each other 
	intersect = [1]

	created = False

	while created != True:

		while sum(intersect) != 0: 

			intersect = []                                                
			if n <= 50: 
				z = [1]*n_nonprimitives + sorted(generate_degree_seq (n_primitives))
			else: 
				z = [1]*n_nonprimitives + sorted(generate_degree_seq2 (n_primitives)) 

			G = nx.configuration_model(z, create_using=nx.Graph)   # configuration model

			while nx.is_connected(G) == False:                     # make sure the graph is connected
				if n <= 50: 
					z = [1]*n_nonprimitives + sorted(generate_degree_seq (n_primitives))
				else: 
					z = [1]*n_nonprimitives + sorted(generate_degree_seq2 (n_primitives))
				
				G = nx.configuration_model(z, create_using=nx.Graph)
			
			nd2 = [n for n, d in G.degree() if d==2]    # number of nodes with degree 2
			# print('d2 nodes', nd2)
			for n in nd2: intersect.append(len(set(list(G.neighbors(n))).intersection(set(nd2)))) 

		# detect self-loops 
		selfloop = list(nx.selfloop_edges(G))
		nd1 = [n for n, d in G.degree() if d==1]        # number of nodes with degree 1

		if selfloop == [] and nd1 == [0,1,2,3]:
			nx.write_edgelist(G, edgelist_filename (settings, sample))
			created = True



##########################################
### run partition                      ### 
##########################################

def run_metis (G, nparts):
	# this function partitions a network into k subgraphs 
	partition = metis.part_graph(G, nparts = nparts, recursive=True)
	return partition

def find_best_n (G, nonprimitives): 
	""" for a given G, obtain the optimal nparts that gives the minimum cuts """
	# remove input and output nodes 
	G_primitive = copy.deepcopy(G)
	for node in nonprimitives:
		G_primitive.remove_node(node)

	min_nodes = 5
	max_nodes = 10
	nparts = list(range(2, int(len(G_primitive.nodes())/min_nodes)+1))

	# find n (block number) that gives the minimum cuts 
	mincuts, part_opt = 100, []
	for n in nparts:
		partition = run_metis (G_primitive, n)
		if partition[0] < mincuts:
			part_opt = partition 
			mincuts = part_opt[0]

	return part_opt

def partition(settings, sample):
	""" 
	partition a graph 
	"""
	inputfile = DAG_edgelist_filename (settings, sample)
	G = nx.Graph()
	G = nx.read_edgelist(inputfile, nodetype=str)

	edge_direct = read_edge_direction (settings, sample)

	nonprimitives = settings[sample]['in_nodes'].split(',') + settings[sample]['out_nodes'].split(',')
	part_opt = find_best_n(G, nonprimitives)

	# write partition solution
	outfile = part_sol_filename (settings, sample)
	G_primitive = copy.deepcopy(G)
	for node in nonprimitives:
		G_primitive.remove_node(node)

	f_out = open(outfile, 'w')
	f_out.write('cut\t'+str(part_opt[0])+'\n')
	for part in range(max(part_opt[1])+1):
		nodeIdx = [a for a, b in enumerate(part_opt[1]) if b == part]
		nodes = [list(G_primitive.nodes())[n] for n in nodeIdx] 
		f_out.write('Partition '+str(part)+'\t')
		f_out.write(','.join([str(n) for n in nodes])+'\n')

	# visualize partition results
	visualize_assignment (G, part_opt, nonprimitives, part_outfig (settings, sample))

	return G, edge_direct, part_opt



##########################################
### assign biological gates            ### 
##########################################

def calc_circuit_score (assignedLib, gateLib, sensorLib):
	logic = load_logic (settings, sample)

	# initialize a score dictionary
	scoreDict = {}
	for node in assignedLib:
		if assignedLib[node]['type'] != 'input':
			scoreDict[node] = {'logic': [], 'output REU': []}

	for i in range(8):
		# print('state ', i)
		visited = 0
		assignedLib_i = copy.deepcopy(assignedLib)   # set up a temporary library

		# first calculate the input REU 
		for node in assignedLib_i:
			if assignedLib_i[node]['type'] == 'input':
				if logic[assignedLib_i[node]['sensor']][i] == 0: 
					assignedLib_i[node]['output REU'] = assignedLib_i[node]['REU OFF']
					assignedLib_i[node]['logic'] = 0
				else: 
					assignedLib_i[node]['output REU'] = assignedLib_i[node]['REU ON']
					assignedLib_i[node]['logic'] = 1
				# print(i, node, assignedLib_i[node]['sensor'], assignedLib_i[node])
				visited += 1

		# calculate the REU of primitive and output node
		for node in assignedLib_i:
			if len(assignedLib_i[node]['in']) == 1:
				assignedLib_i[node]['visited'] = [0]
				assignedLib_i[node]['logic'] = [-1]
			elif len(assignedLib_i[node]['in']) == 2:
				assignedLib_i[node]['visited'] = [0, 0]
				assignedLib_i[node]['logic'] = [-1, -1]
		r = 1

		while visited != len(assignedLib_i.keys()):
			# print('round ###########################################', r)
			for node in assignedLib_i:
				if assignedLib_i[node]['output REU'] == 0:
					# print('node', node)
					# get input REU
					# print('incoming nodes', assignedLib_i[node]['in'])
					in_nodes = assignedLib_i[node]['in']
					for idx, in_node in enumerate(in_nodes):
						# print('input node', in_node, assignedLib_i[in_node]['output REU'], assignedLib_i[in_node]['logic'])
						# if the in node has a calculated output REU 
						if assignedLib_i[in_node]['output REU'] != 0 and assignedLib_i[node]['visited'][idx] != 1:
							# print(in_node, assignedLib_i[in_node]['output REU'])
							assignedLib_i[node]['input REU'].append(assignedLib_i[in_node]['output REU'])
							assignedLib_i[node]['visited'][idx] = 1
							assignedLib_i[node]['logic'][idx] = assignedLib_i[in_node]['logic']

					# print('inputREU', assignedLib_i[node]['input REU'])
					# print(assignedLib_i[node]['visited'])
					# output REU
					if 0 not in assignedLib_i[node]['visited']:
						if assignedLib[node]['type'] != 'output':
							params = assignedLib_i[node]['params']
							x      = sum(assignedLib_i[node]['input REU'])
							# print('x', x)
							# print('params', params)
							assignedLib_i[node]['output REU'] = sigmoid(params[0], params[1], params[2], params[3], x)
							if 1 in assignedLib_i[node]['logic']:
								assignedLib_i[node]['logic'] = 0
							else: 
								assignedLib_i[node]['logic'] = 1
							# print('output REU', assignedLib_i[node]['output REU'])
							# print('logic of gate', node, assignedLib_i[node]['logic'])
							# print('number of gates that have output REU', visited)
						else: 
							assignedLib_i[node]['output REU'] = sum(assignedLib_i[node]['input REU'])
							if 1 in assignedLib_i[node]['logic']:
								assignedLib_i[node]['logic'] = 1
							else: 
								assignedLib_i[node]['logic'] = 0
							# print('done')
							# update score dictionary

						scoreDict[node]['logic'].append(assignedLib_i[node]['logic'])
						scoreDict[node]['output REU'].append(assignedLib_i[node]['output REU'])
						visited += 1
			r += 1

	# calculate score of this permutation
	Smin = 1e3
	for node in assignedLib:
		if assignedLib[node]['type'] == 'output':
			# print(node)
			# print(scoreDict[node]['output REU'])
			# print(scoreDict[node]['logic'])
			maxOFF = max([scoreDict[node]['output REU'][s] for s in range(8) if scoreDict[node]['logic'][s] == 0])
			minON = min([scoreDict[node]['output REU'][s] for s in range(8) if scoreDict[node]['logic'][s] == 1])
			# print('min on', minON, 'max off', maxOFF)
			S = minON/maxOFF
			if S < Smin: 
				Smin = S
	# return the lowest S
	return Smin

def record_library (assignedLib, outfile):
	""" record the library with gate assignments that generate the highest circuit score """
	f_out = open(outfile, 'w')
	f_out.write('node\ttype\tsensor/gate\tparams\tpartition\n')
	for node in assignedLib: 
		if assignedLib[node]['type'] == 'input':
			f_out.write('\t'.join([ node, assignedLib[node]['type'], assignedLib[node]['sensor'], str([assignedLib[node]['REU ON'], assignedLib[node]['REU OFF']]), 'na' ])+'\n')
		elif assignedLib[node]['type'] == 'primitive':
			f_out.write('\t'.join([ node, assignedLib[node]['type'], assignedLib[node]['gate'], str(assignedLib[node]['params']), str(assignedLib[node]['part']) ])+'\n')
		else: 
			f_out.write('\t'.join([ node, assignedLib[node]['type'], 'na', 'na', 'na'])+'\n')

def assign_gates (G, edge_direct, partition):
	# assign biological gates to partitioned graphs
	gateLib = load_gate_info( settings, sample )
	sensorLib = load_sensor_info( settings, sample )

	assignedLib = {}
	# add 'input node' and 'output node' of each gate
	for v in G.nodes():
		assignedLib[v] = {}
		assignedLib[v]['in'], assignedLib[v]['out'] = [], []
		assignedLib[v]['input REU'] = []
		assignedLib[v]['output REU'] = 0
	for e in edge_direct:
		assignedLib[e[0]]['out'].append(e[1])
		assignedLib[e[1]]['in'].append(e[0])	

	# assign input and output nodes
	nonprimitives = settings[sample]['nonprimitives'].split('\t')
	inodes = settings[sample]['in_nodes'].split(',')
	onodes = settings[sample]['out_nodes'].split(',')
	sensors = ['PTac', 'PTet', 'PBAD']
	for idx, v in enumerate(inodes):
		assignedLib[v]['REU ON']  = sensorLib[sensors[idx]]['on']
		assignedLib[v]['REU OFF'] = sensorLib[sensors[idx]]['off']
		assignedLib[v]['sensor']  = sensors[idx]
		assignedLib[v]['type']    = 'input'
	for v in onodes: 
		assignedLib[v]['type']    = 'output'

	# copy another G, remove input and output nodes 
	G_primitive = copy.deepcopy(G)
	for node in inodes+onodes:
		G_primitive.remove_node(node)

	## initialize the gate assignment
	for i in range(0, max(partition[1])+1):
		# print(partition[1])
		nodeIdx = [a for a, b in enumerate(partition[1]) if b == i]
		gnodes = [list(G_primitive.nodes())[n] for n in nodeIdx] 
		# print(gnodes, 'gate nodes within this partition')
		# update partition result
		for v in gnodes:
			assignedLib[v]['part'] = i   # partition of this gate
		# assign repressor gate
		chosengates = random.sample(range(1,13), len(gnodes))
		# print(chosengates)
		for idx, v in enumerate(gnodes): 
			vg = chosengates[idx]
			vg_rbs = [key for key in gateLib if key.split('-')[0] == str(vg)]  # different rbs choices for this repressor gate
			vg = random.choice(vg_rbs) 
			assignedLib[v]['gate'] = vg
			assignedLib[v]['params'] = [float(gateLib[vg]['ymin']), float(gateLib[vg]['ymax']), float(gateLib[vg]['K']), float(gateLib[vg]['n'])]
			assignedLib[v]['type'] = 'primitive'

 	# calculate circuit score of initial assignments
	S0 = calc_circuit_score (assignedLib, gateLib, sensorLib)

	# simulated annealing
	Tmax = 500   # starting temperature
	C = 1e-5

	for i in range(10):  # for 100 trajectories
		SList = []
		Smax = S0
		bestLib = copy.deepcopy(assignedLib)

		for t in range(500000):   # perform simulated annealing
			# print('trial', t)
			# clone a assignedLib
			assignedLib_tmp = copy.deepcopy( assignedLib )
			# print('new tmp assignedLib', assignedLib_tmp)
			# randomly select a node
			g_swap = random.choice(list(G_primitive.nodes()))
			g_part = assignedLib_tmp[g_swap]['part']
			# print(g_swap, g_part, assignedLib_tmp[g_swap]['gate'], 'gate to be swapped')
			# get all available gates (a gate already in this partition, a different gate with a never used repressor from the library)
			# print('gates within the gate library', list(gateLib.keys()))
			used_gates = [assignedLib_tmp[g]['gate'] for g in assignedLib_tmp if assignedLib_tmp[g]['type'] == 'primitive' and assignedLib_tmp[g]['part'] == g_part and g != g_swap]
			# print('other used gate within this partition', used_gates)
			used_repr = [r.split('-')[0] for r in used_gates]
			# print('other used repressor within this partition', used_repr)
			# remove repressors that have already been used by other nodes in this partition
			availgates = list(set(gateLib.keys())-set([g for g in gateLib if g.split('-')[0] in used_repr]))
			availgates.remove(assignedLib_tmp[g_swap]['gate'])
			# add other gates that have been used in this partition
			availgates = availgates + used_gates 
			# print('available gates', availgates)

			# swap two gates g_swap and g_swap_f
			g_swap_t = random.choice(availgates)
			# print('gate swapped to', g_swap_t)
			if g_swap_t in used_gates:  # if swapping two gates within the same partition
				# update the gate and transfer function params of the swapped gate
				g_swap_f = [g for g in assignedLib_tmp if assignedLib_tmp[g]['type'] == 'primitive' and assignedLib_tmp[g]['part'] == g_part and assignedLib_tmp[g]['gate'] == g_swap_t][0]
				assignedLib_tmp[g_swap_f]['gate'] = assignedLib_tmp[g_swap]['gate']
				assignedLib_tmp[g_swap_f]['params'] = [float(gateLib[assignedLib_tmp[g_swap]['gate']]['ymin']), float(gateLib[assignedLib_tmp[g_swap]['gate']]['ymax']), float(gateLib[assignedLib_tmp[g_swap]['gate']]['K']), float(gateLib[assignedLib_tmp[g_swap]['gate']]['n'])]
				# update the transfer function params of the chosen gate
				assignedLib_tmp[g_swap]['gate'] = g_swap_t
				assignedLib_tmp[g_swap]['params'] = [float(gateLib[g_swap_t]['ymin']), float(gateLib[g_swap_t]['ymax']), float(gateLib[g_swap_t]['K']), float(gateLib[g_swap_t]['n'])]
				# print(assignedLib[g_swap])
				# print(assignedLib_tmp[g_swap])
				# print(assignedLib[g_swap_f])
				# print(assignedLib_tmp[g_swap_f])
			else: # if swapping with a gate in the gate library
				assignedLib_tmp[g_swap]['gate'] = g_swap_t
				# update the transfer function params of the chosen gate
				assignedLib_tmp[g_swap]['params'] = [float(gateLib[g_swap_t]['ymin']), float(gateLib[g_swap_t]['ymax']), float(gateLib[g_swap_t]['K']), float(gateLib[g_swap_t]['n'])]
				# print(assignedLib[g_swap])
				# print(assignedLib_tmp[g_swap])

			# calculate the S score, store it 
			S = calc_circuit_score (assignedLib_tmp, gateLib, sensorLib)
			# print('score', S, 'original score', S0)

			# choose to accept or reject the choice
			
			Ti = Tmax*(math.exp(-C*t))
			if Ti > 0:
				try:
					P = math.exp((S-S0)/Ti)
				except OverflowError: 
					P = 2
				# print('temperature', Ti)
				# print('Probability', P)
			else:
				P = math.exp((S-S0)/0.01)

			# append the highest score that occurs to this point
			if S > Smax: 
				Smax = S
				bestLib = copy.deepcopy(assignedLib_tmp)
			SList.append(Smax)    
			# print('highest score that occurs to round', t, SList)

			# if P > 1, accept change 
			if P > 1:
				assignedLib = assignedLib_tmp
				# print('P>1, accepted')
			# if P < 1, accept change based on probability
			else: 
				if random.random() < P:
					assignedLib = assignedLib_tmp
					# print('P<1, accepted')
				# else: 
					# print('P<1, rejected')

			# print('assigned library after round', t, assignedLib)

		# record the max S
		path = settings[sample]['graph_path']+'/n'+settings[sample]['n']+'/DAG'+settings[sample]['n']+'.'+settings[sample]['k']
		if not os.path.exists(path):
			os.mkdir(path)

		try:
			outfile = open(S_trajectory_filename (settings, sample), 'a')
		except FileNotFoundError:
			outfile = open(S_trajectory_filename (settings, sample), 'w')

		outfile.write(str(i)+'\t')
		outfile.write(','.join([str(S) for S in SList])+'\n')
		# record the best library
		record_library (bestLib, opt_gate_assignment_filename (settings, sample, str(i)))



###########################################
# Refine gate assignment results  ## TO BE UPDATED
###########################################

def generate_combinations (n, rlist):
	""" from n choose r elements """
	combs = [list(itertools.combinations(n, r)) for r in rlist]
	return [item for sublist in combs for item in sublist]

def obt_block_feat (G_primitive):
	""" from a partition result, calculate the load in each block """
	block_load, block_gates = {}, {}
	part = [G_primitive.nodes[n]['part'] for n in G_primitive.nodes()]
	for i in range(0, max(part)+1):
		nodes = [a for a, b in enumerate(part) if b == i]
		block_load[i] = sum([G_primitive.nodes[v]['weight'] for v in nodes]) # summ weight in block i 
		block_gates[i] = [G_primitive.nodes[v]['gate'] for v in nodes]       # create a list of all gates in each block
	return block_load, block_gates

def _choose_nodes (part, owblock):
	""" from a partition result, choose 1-3 nodes from each block """
	# from each overweight block, choose 1-3 nodes (combination)
	comboList = []
	for block in owblock: 
		# print('block #', block)
		nodes = [a for a, b in enumerate(part) if b == block] # nodes in this block 
		comboList.append(generate_combinations (nodes, range(1, 4)))     # select n nodes from block to move out
		# print ('chosen nodes from block', comboList)
	# choose one element (combination of nodes) from each list and assemble
	outnodesList = list(itertools.product(*comboList))
	# print ('list of nodes to be moved out ', outnodesList)
	return outnodesList

def create_nodes_dict (part):
	""" assign a dictionary to store the features of inidividual nodes """
	data_dict = {}
	for node, block in enumerate(part):
		data_dict[node] = block
	return data_dict


def create_temp_part (G, nodes_dict, outnodes, blockbs):
	""" replace the block assignment of outnodes in the partition result and check if new partition is OK """
	part_temp = False

	# assign the partition of nodes in the outnodes to new blocks
	for idx, node in enumerate(outnodes): 
		nodes_dict[node] = blockbs[idx]
	# print('newly assigned partition result', nodes_dict.values())

	# obtain features of each block (load, gates)
	block_load, block_gates = obt_block_feat (G, nodes_dict.values())
	# print(block_load)
	block_ow = [i for i in block_load if block_load[i]>12.92]                                               # list of blocks that are overweight
	block_dupgate = [i for i in block_gates.keys() if len(block_gates[i]) != len(set(block_gates[i]))]   # list of blocks that have dupliate gates

	# check if all blocks in partition are within weight limit and have no repetitive gates
	if len(block_ow) == 0 and len(block_dupgate) == 0:
		part_temp = nodes_dict

	return part_temp, block_load


def cal_cut (G, part):
	""" calculate the change in cut size """
	cut = 0
	for edges in G.edges():
		v1, v2 = edges[0], edges[1]
		if part[v1] != part[v2]:       # if v1 and v2 are NOT in the same block
			cut += 1
	return cut


def _move_nodes (G, part_opt, owblock, uwblock, block_load):
	""" for a given block a, move up to 3 nodes to available blocks and compute the gain """
	
	solns, idx = {}, 1
	outList = _choose_nodes (part_opt[1], owblock)      # select nodes from overweight block to move to

	for outnodes in outList:
		outnodes = list(sum(outnodes, ()))
		blockbList = [p for p in itertools.product(uwblock, repeat=len(outnodes))]    # choose target block to move these nodes to
		# for each permutation, update partition results and check if new partition is OK
		# print('outnodes', outnodes)
		# print('new locations/blocks of these nodes', blockbList)
		for blockbs in blockbList: 
			nodes_dict = create_nodes_dict(part_opt[1])
			# print('original partition', nodes_dict.values())
			part_temp, new_load = create_temp_part(G, nodes_dict, outnodes, blockbs)
			if part_temp != False: 
				# if new partition is OK (all blocks have no duplicated gates nor being over weight limit), 
				# calculate the change in cut size 
				cut = cal_cut (G, part_temp)
				# store delta cut size in solutions 
				gates = [G.nodes[v]['gate'] for v in G.nodes()]
				solns[idx] = {'original partition': part_opt[1], 'loads': block_load.values(), 'new partition': list(part_temp.values()), \
								'gates': gates, 'new loads': new_load.values(), 'new cut': cut, 'original cut': part_opt[0]}
				idx += 1
	return solns

def refine_blocks(settings, sample):
	""" after gate assignment, for any block with high toxicity (>12.92), move nodes from that block to others 
	with the goal to minimize cuts"""

	# part_opt = find_best_n (G)
	# G = assign_gates(G, part_opt)   

	G = nx.Graph()
	G = nx.read_edgelist(DAG_edgelist_filename (settings, sample), nodetype=str)
	G = load_gate_assignment (settings, sample, G)

	G_primitive = copy.deepcopy(G)
	nonprimitives = settings[sample]['in_nodes'].split(',') + settings[sample]['out_nodes'].split(',')
	for node in nonprimitives:
		G_primitive.remove_node(node)

	# for the optimal par2tition that yields the minimum cuts, calculate the average weights in each block
	block_load, block_gates = obt_block_feat (G_primitive)
	block_ow = [i for i in block_load if block_load[i]>12.92]  # list of blocks that are overweight
	block_uw = list(set(block_load.keys()) - set(block_ow)) # list of blocks that are underweight
	# check if block_uw can add at least one more nodes
	block_uw = [i for i in block_uw if block_load[i]<11.98]
	gates = [G.nodes[v]['gate'] for v in G.nodes()]
	loads = block_load.values() 

	if len(block_ow) in [2, 3]:
		# print('overloaded blocks', block_ow)
		solns = _move_nodes (G, part_opt, block_ow, block_uw, block_load)
		if solns == {}:
			solns = {-1: {'original partition': part_opt[1], 'original cut': part_opt[0], 'gates': gates, 'loads': loads}}

	else: # no overweight block or more than 3 overweight blocks
		if len(block_ow) == 0:
			solns = {0: {'original partition': part_opt[1], 'original cut': part_opt[0], 'gates': gates, 'loads': loads}}
		else: 
			solns = {-1: {'original partition': part_opt[1], 'original cut': part_opt[0], 'gates': gates, 'loads': loads}}

	# print(solns)
	return G, solns

def record_solns (G, solns, outfile): 
	""" record solutions: idx, original cut, original partition, delta cut, new partitions """
	f_out = open(outfile, 'w')
	f_out.write('Solution\tOri Cut\tOri Load\tOri Partition Result\tGate Assignment\tNew Cut\tNew Load\tNew Partition\n')
	for idx in solns: 
		soln = solns[idx]
		if 'new partition' in soln:
			f_out.write('\t'.join([str(idx), str(soln['original cut']), ','.join([str(i) for i in soln['loads']]), \
									','.join([str(i) for i in soln['original partition']]), ','.join(soln['gates']),\
									str(soln['new cut']), ','.join([str(i) for i in soln['new loads']]), \
									','.join([str(i) for i in soln['new partition']])])+'\n')
		else: 
			f_out.write('\t'.join([str(idx), str(soln['original cut']), ','.join([str(i) for i in soln['loads']]), \
									','.join([str(i) for i in soln['original partition']]), ','.join(soln['gates'])])+'\n')


def compute_refinement_results(k):

	for i in range(2):
		# print('Graph ', i)
		# load G 
		G = nx.read_edgelist('./Random graphs/n'+str(k)+'/G'+str(k)+'.'+str(i)+'.edgelist', nodetype = int)
		G, solns = refine_blocks(G)
		outfile = './Random graphs/n'+str(k)+'/G'+str(k)+'.'+str(i)+'.solutions.txt'
		record_solns (G, solns, outfile)



##########################################
### Plotting modules                   ###
##########################################

def visualize_assignment (G, partition, nonprimitives, outfig): 
	""" visualize the partitioned graph with each color representing a partitioned block """
	pos =nx.spring_layout(G)
	fig = plt.figure(figsize=(5,5))
	ax = fig.add_subplot(1,1,1)
	color = cm.rainbow(np.linspace(0,1,max(partition[1])+1))

	G_primitive = copy.deepcopy(G)
	for node in nonprimitives:
		G_primitive.remove_node(node)

	# draw input and output nodes
	nx.draw_networkx_nodes(G, pos, nodelist = nonprimitives, node_color = 'lightgrey')
	labels = {n:n for n in nonprimitives}
	nx.draw_networkx_labels(G, pos, labels)

	# draw gate nodes
	for i in range(max(partition[1])+1):
		nodeIdx = [a for a, b in enumerate(partition[1]) if b == i]
		nodes = [list(G_primitive.nodes())[n] for n in nodeIdx] 
		print('part', i, nodes)

		nx.draw_networkx_nodes(G, pos, nodelist = nodes, node_color = color[i])
		labels = {n:n for n in nodes}
		nx.draw_networkx_labels(G, pos, labels)

	nx.draw_networkx_edges(G, pos, edgelist=G.edges(), arrowstyle='->')
	ax.axis('off')
	plt.savefig(outfig, dpi=200)
	# plt.show() 

