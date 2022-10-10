# genetic-circuit-partitioning
Partitioning a randomly synthesized graph, or a given network of interest 

First download all dependacies in requirements.txt

Users may provide constraits and file locations in the /runs/settings.txt file
Main functions are written in genetic_partition.py 

Description of wrapper functions in the /bin folder: 

1_generate edgelist.py 
Generate DAG edgelist from Cello's JSON output 

2_graph_partition.py
Perform graph partition with METIS (input graphs: circuit DAGs)

Alternatively, if partitioning a biological network, use 2_graph_partition_bionetworks.py

3_partition_optimize.py
Partition optimization 


