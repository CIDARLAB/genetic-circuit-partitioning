
#	Copyright (C) 2020 by
#	Jing Zhang <jgzhang@bu.edu>, Densmore Lab, Boston University
# 	All rights reserved.
#	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

# Map raw reads

BIN_PATH=/home/unix/jgzhang/genetic-circuit-partitioning/bin

qsub -o ./logs/01_partiton_gate_assignment.py_1_1.out.log -e ./logs/01_partiton_gate_assignment.py_1_1.err.log "python $BIN_PATH/partiton_gate_assignment.py.py -settings ./data/settings.txt -samples 1"
qsub -o ./logs/01_partiton_gate_assignment.py_2_1.out.log -e ./logs/01_partiton_gate_assignment.py_2_1.err.log "python $BIN_PATH/partiton_gate_assignment.py.py -settings ./data/settings.txt -samples 2"
qsub -o ./logs/01_partiton_gate_assignment.py_3_1.out.log -e ./logs/01_partiton_gate_assignment.py_3_1.err.log "python $BIN_PATH/partiton_gate_assignment.py.py -settings ./data/settings.txt -samples 3"
qsub -o ./logs/01_partiton_gate_assignment.py_4_1.out.log -e ./logs/01_partiton_gate_assignment.py_4_1.err.log "python $BIN_PATH/partiton_gate_assignment.py.py -settings ./data/settings.txt -samples 4"
qsub -o ./logs/01_partiton_gate_assignment.py_5_1.out.log -e ./logs/01_partiton_gate_assignment.py_5_1.err.log "python $BIN_PATH/partiton_gate_assignment.py.py -settings ./data/settings.txt -samples 5"
qsub -o ./logs/01_partiton_gate_assignment.py_6_1.out.log -e ./logs/01_partiton_gate_assignment.py_6_1.err.log "python $BIN_PATH/partiton_gate_assignment.py.py -settings ./data/settings.txt -samples 6"
qsub -o ./logs/01_partiton_gate_assignment.py_7_1.out.log -e ./logs/01_partiton_gate_assignment.py_7_1.err.log "python $BIN_PATH/partiton_gate_assignment.py.py -settings ./data/settings.txt -samples 7"
qsub -o ./logs/01_partiton_gate_assignment.py_8_1.out.log -e ./logs/01_partiton_gate_assignment.py_8_1.err.log "python $BIN_PATH/partiton_gate_assignment.py.py -settings ./data/settings.txt -samples 8"
qsub -o ./logs/01_partiton_gate_assignment.py_9_1.out.log -e ./logs/01_partiton_gate_assignment.py_9_1.err.log "python $BIN_PATH/partiton_gate_assignment.py.py -settings ./data/settings.txt -samples 9"