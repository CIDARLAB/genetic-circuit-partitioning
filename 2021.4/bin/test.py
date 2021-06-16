import timeit
from copy import deepcopy
# import genetic_partition_test as gp 
import networkx as nx
import numpy as np
import ujson

d1 = {0:['1', '2', '3'], 1:['2', '5', '9'], 3: ['7', '9', '10']}
t1 = [(0, ['1', '2', '3']), (1, ['2', '5', '9']), (3, ['7', '9', '10'])]

PATH = '/Users/jgzhang/Work/Densmore_lab/genetic-circuit-partitioning/2021.4/'

edgelist = PATH + 'runs/benchmark/4-input-boolean-circuits/276/DAG.edgelist'

# G = nx.read_edgelist (edgelist, nodetype = str, create_using=nx.DiGraph())
# in_nodes, out_nodes, nonprimitives  = gp.get_nonprimitive_nodes(G)

# starttime = timeit.default_timer()
# d2 = ujson.loads(ujson.dumps(d1))
# print("Time exec time is ", timeit.default_timer() - starttime)

# starttime = timeit.default_timer()
# d2 = deepcopy(d1)
# print("Time exec time is ", timeit.default_timer() - starttime)


def convert_tuple_to_dict (t):
	d = {}
	for k, v in t: 
		d[k] = v
	return d

tuple(list(partDict.items()))

d = convert_tuple_to_dict (t1)
