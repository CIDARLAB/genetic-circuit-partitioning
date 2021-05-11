""" this script generates random graphs, caveat of this method is that it is difficult to assign directions after graph is made, 
and does not necessarily resemebles a tree-like structure typically seen in a genetic circuit """


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
from networkx.drawing.nx_agraph import graphviz_layout

###########################################
# Generate random graphs as benchmarks
###########################################

def generate_degree_seq(n):
	# for a given number n, randomly generate n1 number of NOT (degree: 2), 
	# and n2 number of NOR (degree: 3) gates
	seq = [1]
	while sum(seq) %2 != 0:
		seq = [random.randint(2,3) for x in range(n)]
	return seq

def generate_degree_seq2(n):
	# as we know that the ratio of degree 3 nodes and degree 2 nodes is 2.16:1 in order to get a connected graph 
	# we can generate random degree seq satisfying this ratio to get successful graphs
	seq = [1]
	choice_list = [2]*32 + [3]*68
	while sum(seq) %2 != 0:
		seq = random.sample(choice_list, n)
	return seq

def generate_random_graph(n):
	# for a given vertice number n, generate a connected graph 
	n0 = 4                                      # define the number of vertices (input/output) with degree 1 
	if n <= 50: 
		z = [1]*n0 + sorted(generate_degree_seq (n-n0))
	else: 
		z = [1]*n0 + sorted(generate_degree_seq2 (n-n0))
	# G = nx.configuration_model(z)               # configuration model
	G = nx.configuration_model(z, create_using=nx.Graph)
	while nx.is_connected(G) == False:          # make sure the graph is connected
		if n <= 50: 
			z = [1]*n0 + sorted(generate_degree_seq (n-n0))
		else: 
			z = [1]*n0 + sorted(generate_degree_seq2 (n-n0))
		G = nx.configuration_model(z, create_using=nx.Graph)
	return G, Counter(z)


def generate_random_graphs(N, x): 
	# for a given list of vertice number n, generate x number of random graphs
	#fig = plt.figure(figsize=(5,5))
	countList, GList = [], []  
	idx = 1         
	while idx <= 100: 
		print( 'Graph ', idx)
		#ax = fig.add_subplot(1,1,i)
		# generate a random graph in which two d2 nodes are not connected to each other 
		intersect = [1]
		while sum(intersect) != 0: 
			intersect = []
			G, count = generate_random_graph(N)
			nd2 = [n for n, d in G.degree() if d==2]    # number of nodes with degree 2
			#print('d2 nodes', nd2)
			for n in nd2: 
				intersect.append(len(set(list(G.neighbors(n))).intersection(set(nd2)))) 

		# detect self-loops 
		selfloop = list(nx.selfloop_edges(G))
		nd1 = [n for n, d in G.degree() if d==1]        # number of nodes with degree 1

		if selfloop == [] and nd1 == [0,1,2,3]:
			fig = plt.figure()
			nx.write_adjlist(G, './Random graphs/n'+str(N)+'/G'+str(N)+'.'+str(idx+100)+'.adjlist')
			nx.write_edgelist(G, './Random graphs/n'+str(N)+'/G'+str(N)+'.'+str(idx+100)+'.edgelist')
			pos = nx.spring_layout(G)
			nx.draw_networkx_nodes(G, pos, nodelist = G.nodes())
			nx.draw_networkx_edges(G, pos, edgelist = G.edges())
			labels = {n:n for n in G.nodes()}
			nx.draw_networkx_labels(G, pos, labels)
			plt.savefig('./Random graphs/n'+str(N)+'/G'+str(N)+'.'+str(idx+100)+'.pdf')
			idx += 1
		else: 
			print('attempt ', idx, selfloop, nd1)
		# add random weights to vertices according to a normal distribution 
		#for v in G.nodes():
		#	G.node[v]['weight'] = np.random.normal(1, 0.1)
		countList.append(count)
		GList.append(G)
		
		# nodes, weights = zip(*nx.get_node_attributes(G, 'weight').items())
		# nx.draw(G, with_labels = True, node_color=weights, node_size=100, cmap=plt.cm.Purples)
	#plt.show()
	#plt.savefig('n20.png', dpi=300)
	return GList, countList


def node_counter(a, b):
	# for a given number x from a to b, generate random graphs consisting of x nodes with the predefined constraint, 
	# and calculate the average numbers of nodes with degree 2 or 3. 
	f_out = open('degree counts.txt', 'a')
	f_out.write('nodes\tdegree 2\tdegree 3\n')
	for i in range(a, b+1): 
		nd2, nd3 = [], []
		if i <= 70:
			GList, countList = generate_random_graphs(i, 50)
		else: 
			GList, countList = generate_random_graphs(i, 10)
		for count in countList: 
			nd2.append(count[2])
			nd3.append(count[3])
		f_out.write('\t'.join([str(i), str(np.average(nd2)), str(np.average(nd3))])+'\n')
	

def determine_self_loops ():
	fig = plt.figure(figsize=(8,8))
	f_out = open('degree one count.txt', 'w')
	for idx, k in enumerate([20], 1):
		degree1_count = []
		for i in range(1, 101):
			G = nx.read_edgelist('./Random graphs/n'+str(k)+'/G'+str(k)+'.'+str(i)+'.edgelist', nodetype = int)
			degreeList = [G.degree(n) for n in G.nodes()]
			degree1_count.append(degreeList.count(1))
			if degreeList.count(1) != 4:
				print(i)
		f_out.write(str(k)+'\t'+','.join([str(c) for c in degree1_count])+'\n')
		ax = fig.add_subplot(2,2,idx)
		ax.hist(degree1_count)
		ax.set_title(k)
	plt.show()

def check_dag (k, i):
	""" for a given graph, check for acyclic-ness """
	G = nx.read_edgelist('./Random graphs/n'+str(k)+'/DAG'+str(k)+'.'+str(i)+'.edgelist', nodetype = int, create_using=nx.DiGraph())
	return nx.is_directed_acyclic_graph(G)


if __name__ == '__main__': 

	generate_random_graphs(20, 100)
