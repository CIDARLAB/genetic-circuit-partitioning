
#	Copyright (C) 2020 by
#	Jing Zhang <jgzhang@bu.edu>, Densmore Lab, Boston University
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

BIN_PATH='/Users/jgzhang/Work/Densmore_lab/genetic-circuit-partitioning/2021.4/bin'

#python3 $BIN_PATH/3_partition_optimize_test.py -settings ./settings.txt -samples RCA4
 python3 -m cProfile -s tottime $BIN_PATH/3_partition_optimize_test.py -settings ./settings.txt -samples 80727

