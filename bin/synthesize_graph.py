#!/usr/bin/env python

#	Copyright (C) 2020 by
#	Jing Zhang <jgzhang@bu.edu>, Densmore Lab, Boston University
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

# Supporting modules
import argparse
import genetic_partition as gp 

def main():
	# Parse the command line inputs
	parser = argparse.ArgumentParser(description="map_reads")
	parser.add_argument("-settings",  dest="settings",  required=True,  help="settings.txt", metavar="string")
	parser.add_argument("-samples",  dest="samples",  required=True,  help="1,2", metavar="string")
	args = parser.parse_args()
	# Run the command
	samples = args.samples.split(',')
	settings = gp.load_settings(args.settings)
	for s in samples:
		generate_random_graph (settings, s)

if __name__ == "__main__":
	main()