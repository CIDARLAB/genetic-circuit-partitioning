
#	Copyright (C) 2020 by
#	Jing Zhang <jgzhang@bu.edu>, Densmore Lab, Boston University
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

BIN_PATH='/home/ubuntu/genetic-circuit-partitioning/2021.4/bin'

python3 $BIN_PATH/3_partition_optimize_bionetworks.py -settings ./settings_bionetwork.txt -samples SWG_n50_k3_p0.5,SWG_n50_k4,SWG_n50_k4_p0.5,SWG_n50

