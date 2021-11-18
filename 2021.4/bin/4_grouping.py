"""
This script aims to group 
"""
import argparse
import genetic_circuit_partition as gp 
from pycallgraph import PyCallGraph
from pycallgraph.output import GraphvizOutput
import os
import datetime

def main():

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

		# load graph 
		G   = gp.load_graph_undirected (settings, s)
		DAG = gp.load_graph (settings, s)

		gp.order_BFS (DAG, primitive_only, S_bounds)


if __name__ == "__main__":
	main()



