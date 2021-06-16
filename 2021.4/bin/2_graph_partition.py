#!/usr/bin/env python

#	Copyright (C) 2021 by
#	Jing Zhang <jgzhang@bu.edu>, Densmore Lab, Boston University
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

# Supporting modules
import argparse
import genetic_partition as gp 
import os
import shutil

def main():

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

		# obtain user-defined params (TO BE MOVED TO MAIN)
		S_bounds        = settings[s]['S_bounds'].split(',')
		target_n        = settings[s]['target_n'].split(',')
		primitive_only  = settings[s]['primitive_only']
		high_constraint = settings[s]['high_constraint']
		low_constraint  = settings[s]['low_constraint']
		out_path        = settings[s]['output_path']

		# load graph 
		G   = gp.load_graph_undirected (settings, s)
		DAG = gp.load_graph (settings, s)

		# if primitive only is true, input and output nodes will be removed from list of vertices to be partitioned
		if primitive_only:
			in_nodes, out_nodes, nonprimitives  = gp.get_nonprimitive_nodes (DAG)
			G_primitive = gp.get_G_primitive (G, nonprimitives)
		else:
			G_primitive   = copy.deepcopy(G)
			nonprimitives = []

		# make output directories
		if os.path.exists(out_path):
			shutil.rmtree(out_path)
			os.mkdir(out_path)
		else:
			os.mkdir(out_path)

		# partition the graph into n parts 
		min_nodes = int (S_bounds[0])
		max_nodes = int (S_bounds[1])

		if target_n == ['']:  # if user doesn't supply a target number that the graph should be partitioned into 
			nparts = list(range(int(len(G_primitive.nodes())/max_nodes)+1, int(len(G_primitive.nodes())/min_nodes)+1))
			os.mkdir (out_path + '/nparts')
			for n in nparts: 
				outdir   = out_path + '/nparts/' + str(n)
				os.mkdir (outdir)
				part_opt = gp.partition_nparts_wrapper (G_primitive, n, outdir)
				gp.visualize_assignment_graphviz (DAG, part_opt, nonprimitives, primitive_only, outdir, 0)

		else:                 # if user provides a target number
			os.mkdir(out_path+'/nparts/')
			os.mkdir(out_path+'/nparts/'+target_n)
			part_opt = gp.partition_nparts_wrapper (G_primitive, int(target_n), out_path+'/nparts/'+target_n)
			gp.visualize_assignment_graphviz (DAG, part_opt, nonprimitives, primitive_only, out_path+'/nparts/'+target_n, 0)


if __name__ == "__main__":
	main()


	