import genetic_partition_test as gp 
import os
import networkx as nx
import random
import itertools

PATH = '/Users/jgzhang/Work/Densmore_lab/Partition/code_version/v2/genetic-circuit-partitioning/2021.4/runs/benchmark/'

def load_graph (edgelist):
	G = nx.read_edgelist (edgelist, nodetype = str, create_using=nx.DiGraph())
	return G

Ninput = [4,5,6,7,8]

node_nums = []

for n in Ninput: 
	bm_dir = PATH + str(n) + '-input-boolean-circuits/'
	for s in os.listdir(bm_dir):
		if s.isdigit():
			bm = bm_dir + s + '/DAG.edgelist'
			if os.path.exists(bm):
				G = load_graph (bm)
				in_nodes, out_nodes, nonprimitives  = gp.get_nonprimitive_nodes (G)
				node_nums.append((n, s, len(list(G.nodes()))-n))

# print(node_nums)

# for each range (10~110), randomly choose 5 circuits  
samples = []

for i in range(10, 110, 10):
	nrange = [i, i+9]
	samples_to_choose = [s for s in node_nums if s[-1]>=nrange[0] and s[-1] <= nrange[1]]
	# print(i, samples_to_choose)
	samples.append(random.sample(samples_to_choose, 5))
	print(i, samples)


samples = list(itertools.chain(*samples))
f_out = open('settings_boolean_circuits.txt', 'w')
for s in samples: 
	ninput = str(s[0])
	name = s[1]
	print(ninput)
	print(name)
	graph = ninput + 'input_'+name
	graph_path = './benchmark/'+ ninput +'-input-boolean-circuits/'+ name
	lib_path = './lib'
	primitive_only = 'True'
	S_bounds = '1,10'
	target_n = ''
	maxNodes = '3'
	hc = '2,1'
	lc = '3'
	loop_free = 'TRUE'
	
	t  = '3'
	output_path = './results/'+ ninput +'-input-boolean-circuits/'+ name
	f_out.write('\t'.join([graph, graph_path, lib_path, primitive_only, S_bounds, target_n, maxNodes, hc, lc, t, output_path])+'\n')




