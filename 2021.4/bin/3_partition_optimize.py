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
		print (settings[s])
		# obtain user-defined params 
		S_bounds        = settings[s]['S_bounds'].split(',')
		target_n        = settings[s]['target_n'].split(',')
		primitive_only  = settings[s]['primitive_only']
		maxNodes        = int(settings[s]['max_nodes_to_shift'])
		high_constraint = settings[s]['high_constraint'].split(',')
		low_constraint  = settings[s]['low_constraint'].split(',')
		loop_free       = settings[s]['loop_free']
		priority        = settings[s]['priority']
		trajectories    = int(settings[s]['trajectories'])
		out_path        = settings[s]['output_path']
		timestep = 100000

		# load graph 
		G   = gp.load_graph_undirected (settings, s)
		DAG = gp.load_graph (settings, s)

		outdir = out_path + '/nparts/'
		order = gp.rank_connectivity (DAG, primitive_only, outdir)
		print('median degree of connectivity in subgraphs', order)

		order = gp.rank_constraint_met (DAG, primitive_only, high_constraint, loop_free, outdir)
		print('percent of cells with unmet constraint under high constraint', order)

		order = gp.rank_constraint_met (DAG, primitive_only, low_constraint, loop_free, outdir)
		print('percent of cells with unmet constraint under low constraint', order)

		# optimize graph partition results to satisfy constraints
		if target_n == ['']:
			nparts = [int(n) for n in os.listdir(outdir)]
			if len(nparts) >= 11: 
				nparts = sorted(nparts)[0:10]
			for npart in sorted(nparts):
				print('target npart', npart)
				part_sol = outdir + str(npart) + '/part_solns.txt'
				cut, partDict = gp.load_metis_part_sol (part_sol)
				if len(list(G.nodes()))>=20:
					print('optimizing by subnetworks')
					print('optimizing with low constraint')
					gp.optimize_signal_subnetwork (DAG, primitive_only, S_bounds, cut, partDict, maxNodes, 3, True, low_constraint, loop_free, priority, timestep, trajectories, outdir + str(npart) +'/optimized_lc/')
					print('execution time', datetime.datetime.now()-begin_time)
					# print('optimizing with high constraint')
					# gp.optimize_signal_subnetwork (DAG, primitive_only, S_bounds, cut, partDict, maxNodes, 5, high_constraint, loop_free, priority, trajectories, outdir + str(npart) +'/optimized_hc/')
					# print('execution time', datetime.datetime.now()-begin_time)
				else:                                
					print('brute-force optimizing')
					print('optimizing with low constraint')
					gp.optimize_signal_bruteforce (DAG, primitive_only, S_bounds, cut, partDict, maxNodes, low_constraint, outdir + str(npart) +'/optimized_lc/')
					print('execution time', datetime.datetime.now()-begin_time)
					# print('optimizing with high constraint')
					# gp.optimize_signal_bruteforce (DAG, primitive_only, S_bounds, cut, partDict, maxNodes, high_constraint, outdir + str(npart) +'/optimized_hc/')
					# print('execution time', datetime.datetime.now()-begin_time)
			# gp.determine_best_solution (DAG, primitive_only, high_constraint, low_constraint, outdir)

		else:
			for n in target_n:
				print('target npart', n)
				outdir = out_path + '/nparts/'
				part_sol = outdir + n + '/part_solns.txt'
				cut, partDict = gp.load_metis_part_sol (part_sol)
				if len(list(G.nodes()))>=25:
					print('optimizing by subnetworks')
					print('optimizing with low constraint')
					gp.optimize_signal_subnetwork (DAG, primitive_only, S_bounds, cut, partDict, maxNodes, 3, True, low_constraint, loop_free, priority, timestep, trajectories, outdir + n + '/optimized_lc/')
					# print('optimizing with high constraint')
					# gp.optimize_signal_subnetwork_qs (DAG, primitive_only, S_bounds, cut, partDict, maxNodes, high_constraint, loop_free, priority, trajectories, outdir + n + '/optimized_hc/')
				else:
					print('brute-force optimizing')
					gp.optimize_signal_bruteforce (DAG, primitive_only, S_bounds, cut, partDict, maxNodes, low_constraint, outdir + n + '/optimized_lc/')
					gp.optimize_signal_bruteforce (DAG, primitive_only, S_bounds, cut, partDict, maxNodes, low_constraint, outdir + n + '/optimized_hc/')
			
				# # if optimizing based on previous results
				# iteration = 2
				# solDict = gp.load_opt_part_sol (outdir + n + '/optimized_lc/part_solns.txt')
				# T = int(solDict[iteration]['T'])
				# cut = int(solDict[iteration]['cut'])
				# part = solDict[iteration]['part']
				# gp.optimize_signal_subnetwork (DAG, primitive_only, S_bounds, cut, part, maxNodes, low_constraint, loop_free, priority, trajectories, outdir + n + '/optimized_lc/'+str(iteration)+'/')

				## if partitioning into target_n doesn't have af solution, try splitting cells 
				# gp.split_cells (DAG, primitive_only, S_bounds, cut, partDict, maxNodes, low_constraint, loop_free, priority, trajectories, outdir + n + '/optimized_lc/')

			# gp.determine_best_solution (DAG, primitive_only, high_constraint, low_constraint, outdir)

		


	print('TOTAL execution time', datetime.datetime.now()-begin_time)

if __name__ == "__main__":
	main()


	