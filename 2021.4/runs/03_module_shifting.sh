
#	Copyright (C) 2020 by
#	Jing Zhang <jgzhang@bu.edu>, Densmore Lab, Boston University
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

BIN_PATH='/Users/jgzhang/Work/Densmore_lab/Partition/code_version/v2/genetic-circuit-partitioning/2021.4/bin'

python3 $BIN_PATH/3_partition_optimize.py -settings ./settings.txt -samples 6input_15,6input_30,6input_26,6input_72,6input_5,6input_67,6input_14,7input_59,8input_46,7input_87,6input_54,7input_36,7input_26,8input_80,7input_65,7input_1,8input_6,7input_80,7input_95,7input_83,7input_72,7input_71,7input_33,7input_15,7input_45,7input_69,7input_54,7input_89,7input_90,8input_70,7input_19,7input_2
# python3 -m cProfile -s tottime $BIN_PATH/3_partition_optimize_test.py -settings ./settings.txt -samples 80727

