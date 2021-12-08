
#	Copyright (C) 2020 by
#	Jing Zhang <jgzhang@bu.edu>, Densmore Lab, Boston University
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

BIN_PATH='/home/ubuntu/genetic-circuit-partitioning/2021.4/bin'

python3 $BIN_PATH/3_partition_optimize_bionetworks.py -settings ./settings_bionetwork.txt -samples SFG_n20,SFG_n20_0.35_0.1_0.55,SFG_n20_0.55_0.1_0.35,SFG_n20_0.65_0.1_0.25,SFG_n20_0.75_0.1_0.15

