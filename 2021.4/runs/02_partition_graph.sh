
#	Copyright (C) 2020 by
#	Jing Zhang <jgzhang@bu.edu>, Densmore Lab, Boston University
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

# Perform metis graph partition

BIN_PATH='/home/ubuntu/genetic-circuit-partitioning/2021.4/bin'

python3 $BIN_PATH/2_graph_partition.py -settings ./settings_bionetwork.txt -samples SFG_n90_0.3_0.05_0.65_r4,SFG_n90_0.3_0.05_0.65,SFG_n100_0.35_0.1_0.55,SFG_n100_0.35_0.1_0.55_r2,SFG_n100_0.35_0.1_0.55_r3,SFG_n100_0.35_0.1_0.55_r4