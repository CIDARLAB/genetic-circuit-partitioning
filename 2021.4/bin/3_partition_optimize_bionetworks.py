#!/usr/bin/env python

#	Copyright (C) 2021 by
#	Jing Zhang <jgzhang@bu.edu>, Densmore Lab, Boston University
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

# Supporting modules
import argparse
import genetic_partition_test as gp 
from pycallgraph import PyCallGraph
from pycallgraph.output import GraphvizOutput
import os
import datetime
import copy

def main():

	graphviz = GraphvizOutput()
	graphviz.output_file = 'partition_RCA4.png'

	begin_time = datetime.datetime.now()

	# Parse the command line inputs
	parser = argparse.ArgumentParser(description="perform graph partition using metis")
	parser.add_argument("-settings",  dest="settings",  required=True,  help="settings.txt", metavar="string")
	parser.add_argument("-samples",  dest="samples",  required=True,  help="1,2", metavar="string")
	args = parser.parse_args()
	
	# Run the command
	samples = args.samples.split(',')
	settings = gp.load_settings(args.settings)
	for s in samples:
		print ('Processing sample', s)
		# obtain user-defined params 
		S_bounds        = settings[s]['S_bounds'].split(',')
		primitive_only  = settings[s]['primitive_only']
		maxNodes        = int(settings[s]['max_nodes_to_shift'])
		priority        = settings[s]['priority']
		loop_free       = 'FALSE'
		trajectories    = int(settings[s]['trajectories'])
		out_path        = settings[s]['output_path']

		# load graph 
		G   = gp.load_graph_undirected (settings, s)
		DAG = gp.load_graph (settings, s)

		outdir = out_path + '/nparts/'
		connectivity_order = gp.rank_connectivity (DAG, primitive_only, outdir)
		lowest_connectivity = min([c for (n, c) in connectivity_order])
		print('median degree of connectivity in subgraphs', connectivity_order)
		print('lowest median degree of connectivity in subgraph', lowest_connectivity)

		# if primitive only is true, input and output nodes will be removed from list of vertices to be partitioned
		if primitive_only == 'TRUE':
			in_nodes, out_nodes, nonprimitives  = gp.get_nonprimitive_nodes (DAG)
			G_primitive = gp.get_G_primitive (G, nonprimitives)
		elif primitive_only == 'FALSE':
			G_primitive   = copy.deepcopy(G)
			nonprimitives = []

		# optimize graph partition results to satisfy constraints
		nparts = list(range(int(len(G_primitive.nodes())/int(S_bounds[1]))+1, int(len(G_primitive.nodes())/int(S_bounds[0]))+1))
		for npart in sorted(nparts):
			print('target n', npart)
			median_connectivity = [c for (n, c) in connectivity_order if n == npart]
			print('median connectivity', median_connectivity)
			part_sol = outdir + str(npart) + '/part_solns.txt'
			cut, partDict = gp.load_metis_part_sol (part_sol)
			# gp.optimize_signal_subnetwork (DAG, primitive_only, S_bounds, cut, partDict, maxNodes, str(lowest_connectivity), loop_free, priority, trajectories, outdir + str(npart) + '/optimized/')



			# for npart in os.listdir(outdir):
			# 	if npart.isdigit():
			# 		print('target npart', npart)
			# 		part_sol = outdir + npart + '/part_solns.txt'
			# 		cut, partDict = gp.load_metis_part_sol (part_sol)
			# 		if len(list(G.nodes()))>=15:
			# 			print('optimizing by subnetworks')
			# 			print('optimizing with low constraint')
			# 			gp.optimize_signal_subnetwork (DAG, primitive_only, S_bounds, cut, partDict, maxNodes, low_constraint, loop_free, priority, trajectories, outdir + npart +'/optimized_lc/')
			# 			print('execution time', datetime.datetime.now()-begin_time)
			# 			# print('optimizing with high constraint')
			# 			# gp.optimize_signal_subnetwork (DAG, primitive_only, S_bounds, cut, partDict, maxNodes, high_constraint, loop_free, priority, trajectories, outdir + npart +'/optimized_hc/')
			# 			# print('execution time', datetime.datetime.now()-begin_time)
			# 		else:                                
			# 			print('brute-force optimizing')
			# 			print('optimizing with low constraint')
			# 			gp.optimize_signal_bruteforce (DAG, primitive_only, S_bounds, cut, partDict, maxNodes, low_constraint, outdir + npart +'/optimized_lc/')
			# 			print('execution time', datetime.datetime.now()-begin_time)
			# 			print('optimizing with high constraint')
			# 			gp.optimize_signal_bruteforce (DAG, primitive_only, S_bounds, cut, partDict, maxNodes, high_constraint, outdir + npart +'/optimized_hc/')
			# 			print('execution time', datetime.datetime.now()-begin_time)
			# gp.determine_best_solution (DAG, primitive_only, high_constraint, low_constraint, outdir)


	print('TOTAL execution time', datetime.datetime.now()-begin_time)

if __name__ == "__main__":
	main()


	