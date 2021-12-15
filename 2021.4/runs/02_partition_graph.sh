
#	Copyright (C) 2020 by
#	Jing Zhang <jgzhang@bu.edu>, Densmore Lab, Boston University
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

# Perform metis graph partition

BIN_PATH='/home/ubuntu/genetic-circuit-partitioning/2021.4/bin'

python3 $BIN_PATH/2_graph_partition_bionetworks.py -settings ./settings_bionetwork.txt -samples SWG_n60_k4_p0.3,SWG_n60_k4_p0.4,SWG_n70_k4_p0.3,SWG_n70_k4_p0.4,SWG_n80_k4_p0.3,SWG_n80_k4_p0.4,SWG_n90_k4_p0.3,SWG_n90_k4_p0.4,SWG_n100_k4_p0.3,SWG_n100_k4_p0.4