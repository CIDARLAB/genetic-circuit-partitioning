
#	Copyright (C) 2020 by
#	Jing Zhang <jgzhang@bu.edu>, Densmore Lab, Boston University
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

# Perform metis graph partition

BIN_PATH='/home/ubuntu/genetic-circuit-partitioning/2021.4/bin'

python3 $BIN_PATH/2_graph_partition_bionetworks.py -settings ./settings_bionetwork.txt -samples RG_n100_p0.01,RG_n100_p0.011,RG_n100_p0.012,RG_n100_p0.013,RG_n100_p0.014