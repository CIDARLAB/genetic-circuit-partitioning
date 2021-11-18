
import genetic_circuit_partition as gp 
import os
import metis

path = '/Users/jgzhang/Work/Densmore_lab/Partition/code_version/v2/genetic-circuit-partitioning/2021.4/'

settings = gp.load_settings (path + 'runs/settings.txt')

s = 'md5Core'

G   = gp.load_graph_undirected (settings, s)
DAG = gp.load_graph (settings, s)

in_nodes, out_nodes, nonprimitives  = gp.get_nonprimitive_nodes (DAG)
G_primitive = gp.get_G_primitive (G, nonprimitives)

nparts = list(range(30,47))
for n in nparts: 
	outdir   = path + 'runs/results/electronic-circuits/md5Core/nparts/' + str(n)
	if os.path.exists(outdir) == False:
		os.mkdir(outdir)
	# part_opt = gp.partition_nparts_wrapper (G_primitive, n, outdir)
	part_opt = metis.part_graph(G, nparts = n, recursive=True)
	outfile = outdir +'/part_solns.txt'
	f_out = open(outfile, 'w')
	f_out.write('cut\t'+str(part_opt[0])+'\n')
                                                                                                                              
	for part in range(max(part_opt[1])+1):
		nodeIdx = [a for a, b in enumerate(part_opt[1]) if b == part]
		nodes = [list(G.nodes())[node] for node in nodeIdx] 
		f_out.write('Partition '+str(part)+'\t')
		f_out.write(','.join([str(node) for node in nodes])+'\n')
# 	gp.visualize_assignment_graphviz (DAG, part_opt, nonprimitives, primitive_only, outdir, 0)