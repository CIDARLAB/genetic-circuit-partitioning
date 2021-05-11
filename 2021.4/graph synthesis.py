import re
import networkx as nx
import matplotlib.pyplot as plt
import os
from networkx.drawing.nx_agraph import graphviz_layout
import copy

def read_json(inputfile):
	lines = [open(inputfile, 'r').read().strip("\n")][0].split('\n')
	ports, gates = {}, {}
	for idx, line in enumerate(lines):
		line = line.strip()
		if line.startswith('"ports"'):
			p_s = idx
			searchlines = lines[idx+1:]
			for i, sl in enumerate(searchlines, idx):
				if sl.strip().startswith('"cells"'):
					p_e = i+1
		if line.startswith('"cells"'):
			g_s = idx
			searchlines = lines[idx+1:]
			for i, sl in enumerate(searchlines, idx):
				if sl.strip().startswith('"netnames"'):
					g_e = i
	# get information of inputs and outputs
	spacer = [idx+p_s+1 for idx, line in enumerate(lines[p_s+1:p_e]) if ': {' in line.strip()]
	for i, v in enumerate(spacer):
		# get names
		s = lines[v].strip()
		s = re.search('"(.*)"', s)
		el = s.group(1)
		ports[el] = {}
		# get directions 
		s = lines[v+1].strip()
		s = re.search('"direction": "(.*)"', s)
		direction = s.group(1)
		ports[el]['direction'] = direction
		# get bits
		s = lines[v+2].strip()
		bits = s.split('[')[1].split(']')[0].strip()
		ports[el]['bits'] = int(bits)
	# get information of gates 
	spacer = [idx+g_s+1 for idx, line in enumerate(lines[g_s+1:g_e]) if '$abc$' in line.strip()]
	for i, v in enumerate(spacer): 
		# get names
		s = int(lines[v].strip().split('"')[1].split('$')[-1])
		gates[s] = {}
		gates[s]['input'] = {}
		gates[s]['output'] = {}
		# search for attributes of this gate
		if i != len(spacer)-1:
			searchlines = lines[v:spacer[i+1]]
		else: 
			searchlines = lines[v:]
		for sl in searchlines:
			# get gate type
			if sl.strip().startswith('"type"'):
				gatetype = re.search('_(.*)_', sl.strip())
				gates[s]['type'] = gatetype.group(1)
			# get input(s)
			if sl.strip().startswith('"A": [') or sl.strip().startswith('"B": ['):
				port = re.search('"(.*)"', sl).group(1)
				bits = sl.split('[')[1].split(']')[0].strip()
				gates[s]['input'][port] = int(bits)
			# get output 
			if sl.strip().startswith('"Y": ['):
				port = re.search('"(.*)"', sl).group(1)
				bits = sl.split('[')[1].split(']')[0].strip()
				gates[s]['output'][port] = int(bits)
	return ports, gates

def synthesize_graph (ports, gates, outdir, t):
	G = nx.DiGraph()
	# start from the output, add edges
	edges = []

	for p in ports:
		if ports[p]['direction'] == 'output':
			b = ports[p]['bits']
			for g in gates:
				if b == gates[g]['output']['Y']:
					edges.append((g, p))
					# print('output', (g,p))

	for p in ports:
		if ports[p]['direction'] == 'input':
			b = ports[p]['bits']
			for g in gates:
				if b == gates[g]['input']['A']:
					edges.append((p, g))
					# print('input', (p,g))
				if gates[g]['type'] == 'NOR':
					if b == gates[g]['input']['B']:
						edges.append((p, g))
						# print('input', (p,g))

	for g in gates:
		op = gates[g]['output']['Y']
		for sg in gates: 
			if gates[sg]['type'] == 'NOT':
				gin = [gates[sg]['input']['A']]
			else:
				gin = [gates[sg]['input']['A'], gates[sg]['input']['B']]
			if op in gin:
				edges.append((g, sg))
				# print('internal', (g, sg))

	for e in edges:
		G.add_edge(*e)


	# npnodes = ['a', 'b', 'c', 'd']
	# G_primitive = copy.deepcopy(G)
	# for node in npnodes:
	# 	try:
	# 		G_primitive.remove_node(node)
	# 	except nx.exception.NetworkXError:
	# 		pass


	# write graph
	nx.write_adjlist(G, outdir+'/DAG.adjlist')
	nx.write_edgelist(G, outdir+'/DAG.edgelist')

	# plot graph
	fig = plt.figure(figsize=(5,5))
	ax = fig.add_subplot(1,1,1)
	nx.nx_agraph.write_dot(G, outdir+'/DAG.dot')
	pos = graphviz_layout(G, prog='dot')
	# nx.draw(G, pos, with_labels=True, node_size=1, font_size=0.5, arrowsize=1, width=0.5, node_color='#b18ea6')
	nx.draw(G, pos, with_labels=True, node_color='#b18ea6')
	plt.savefig(outdir+'/DAG.pdf', dpi=200)
	plt.show()

if __name__ == '__main__': 

	for i in ['374']:
		path = './examples/'
		file = path+i+'/4input_'+i+'.json'
		dag_file = path+i+'/DAG.edgelist'
		# if os.path.exists(file) and not os.path.exists(dag_file):
		if os.path.exists(file):
			try:
				print(i)
				ports, gates = read_json(file)
				synthesize_graph (ports, gates, path+i, t=10)
			except (UnboundLocalError, ValueError) as e: 
				print('error in file', i)
