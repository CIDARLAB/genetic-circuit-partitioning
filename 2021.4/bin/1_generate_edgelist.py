#!/usr/bin/env python

#	Copyright (C) 2021 by
#	Jing Zhang <jgzhang@bu.edu>, Densmore Lab, Boston University
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

# Supporting modules
import argparse
import genetic_partition as gp 

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
		print('sample ', s)
		graph_path = settings[s]['graph_path']
		JSON_file  = graph_path + '/' + s + '.json'
		# load json file 
		ports, gates = gp.read_json(JSON_file)
		gp.synthesize_graph (ports, gates, graph_path, t=10)
		

if __name__ == "__main__":
	main()


	