import numpy as np
import networkx as nx
import copy 
import matplotlib.pyplot as plt
from networkx.drawing.nx_agraph import graphviz_layout
import seaborn as sns
import random
import os


##########################################################
# Construct a hexagonal grid
##########################################################

def initialize_hexmap (n):
	""" initialize a graph consisting of n*n/2 nodes and connect them in hexagon layout, and associate the pos as coordinates """
	hexmap = nx.Graph()
	# add nodes 
	nodelist = np.arange(1, n**2/2 + 1)
	# associate (x, y) coordinate 
	for iy, y in enumerate(range(n)): # for each row
		for x in range(n): # for each col 
			if x%2 != 0 and y%4 in [0,3]:
				node = int( iy * n/2 + (x//2+1) )
				hexmap.add_node(node, coor=[x,y])
			if x%2 == 0 and y%4 in [1,2]:
				node = int( iy * n/2 + (x//2+1) )
				hexmap.add_node(node, coor=[x,y])

	# add edge
	for i in hexmap.nodes():
		icoor = hexmap.nodes[i]['coor']
		if (icoor[0]%2 == 1 and icoor[1]%4 == 0) or (icoor[0]%2 == 0 and icoor[1]%4 == 2):
			nb_coor = [[icoor[0]-1, icoor[1]+1], [icoor[0]+1, icoor[1]+1], [icoor[0], icoor[1]-1]]
			for j in hexmap.nodes():
				jcoor = hexmap.nodes[j]['coor']
				if jcoor in nb_coor: 
					if (i, j) not in hexmap.edges(): hexmap.add_edge(i, j)
		elif (icoor[0]%2 == 1 and icoor[1]%4 == 3) or (icoor[0]%2 == 0 and icoor[1]%4 == 1):
			nb_coor = [[icoor[0]-1, icoor[1]-1], [icoor[0]+1, icoor[1]-1], [icoor[0], icoor[1]+1]]
			for j in hexmap.nodes():
				jcoor = hexmap.nodes[j]['coor']
				if jcoor in nb_coor:
					if (i, j) not in hexmap.edges(): hexmap.add_edge(i, j)

	return hexmap

def assign_hexmap_color (hexmap, n): 
	""" assign colors to each node on the hexagonal grid, such that each node is surrounded 
	by 3 differently colored nodes """
	for y in range(n): # for each row
		for x in range(n):
			if y%4 == 0:
				if x%4 == 1: assign = 4
				elif x%4 == 3: assign = 1
			elif y%4 == 1:
				if x%4 == 0: assign = 3
				elif x%4 == 2: assign = 2
			elif y%4 == 2:
				if x%4 == 0: assign = 2
				elif x%4 == 2: assign = 3
			elif y%4 == 3:
				if x%4 == 1: assign = 1
				elif x%4 == 3: assign = 4
			for i in hexmap.nodes():
				icoor = hexmap.nodes[i]['coor']
				if [x, y] == icoor: hexmap.nodes[i]['colorcode'] = assign

def visualize_hexmap(hexmap, outfile):
	""" visualize the color assignment on the hexagonal map """

	# create color dictionary
	color = sns.color_palette('hls', n_colors=4)
	colordict = {}
	colorcode = np.arange(1, 5)
	for k in colorcode:
		colordict[k] = color[k-1]

	# draw nodes
	X, Y = [], []
	for node in hexmap.nodes():
		coor = hexmap.nodes[node]['coor']
		plt.plot([coor[0]], [coor[1]], markersize=2, marker='o', color=colordict[hexmap.nodes[node]['colorcode']])
		plt.text(coor[0], coor[1], node)

	# draw edges
	for edge in hexmap.edges():
		n0, n1 = edge[0], edge[1]
		c0, c1 = hexmap.nodes[n0]['coor'], hexmap.nodes[n1]['coor']
		xlist = [c0[0], c1[0]]
		ylist = [c0[1], c1[1]]
		plt.plot(xlist, ylist, 'k', zorder=0)

	plt.box(on=None)
	plt.savefig(outfile+'.pdf', dpi=200)
	plt.show()


##########################################################
# Generate or load a graph 
##########################################################

def generate_graph (n, outfile):
	""" generate a tree """
	G = nx.Graph()
	G.add_node(1) # initialize root
	for i in range(n):
		# randomly sample from nodes that has less than 3 neighbors
		parent = 0
		while parent == 0:
			parent_candidate = random.choice(list(G.nodes()))
			if len(list(G.neighbors(parent_candidate))) < 3:
				parent = parent_candidate
		# print('selecting parent node', parent)
		avail_num_children = 3 - len(list(G.neighbors(parent_candidate)))
		num_children = random.choice(np.arange(1, avail_num_children+1))
		# print('available number for children', avail_num_children)
		# print('randomly generate', num_children, 'children')
		for j in range(num_children):
			child = max(list(G.nodes()))+1
			G.add_node(child)
			G.add_edge(parent, child)
			# print('add edge between', parent, child)
	# write edgelist 
	nx.write_edgelist(G, outfile + '.edgelist')
	return G

def load_graph (dag_file):
	""" load a DIRECTED graph """
	G = nx.read_edgelist(dag_file, nodetype = str, create_using=nx.DiGraph())
	innodes, outnodes = [], []
	for node in G.nodes():
		indegree = G.in_degree(node)
		outdegree = G.out_degree(node)
		if outdegree == 0:
			outnodes.append(node)
		if indegree == 0:
			innodes.append(node)
	return G, innodes, outnodes

def load_undirected_graph (dag_file):
	""" load an UNDIRECTED graph """
	G = nx.read_edgelist(dag_file, nodetype = int)
	return G


def return_prim_graph (dag, innodes, outnodes):
	""" remove input and output nodes from a graph """
	G_primitive = copy.deepcopy(dag)
	nonprims = innodes+outnodes
	for node in nonprims:
		if node in dag.nodes():
			G_primitive.remove_node(node)
	return G_primitive


def visualize_graph_graphviz (G, outfile): 
	""" visualize a network """
	fig = plt.figure(figsize=(4,5))
	ax = fig.add_subplot(1,1,1)

	nx.nx_agraph.write_dot(G, outfile+'.dot')
	pos = graphviz_layout(G, prog='dot')
	nx.draw(G, pos, nodelist = G.nodes(), with_labels=True, node_color='lightgrey')

	plt.savefig(outfile+'.pdf', dpi=200)
	plt.show()


##########################################################
# Assign color to each cell in the graph
##########################################################

def Diff(li1, li2):
	""" get the difference between two lists """
	list_dif = [i for i in li1 + li2 if i not in li1 or i not in li2]
	return list_dif


def assign_colors (G, source):
	""" use breath first search to assign colors """
	colorcode = np.arange(1, 5)
	# start from a node that has no outgoing edge
	if source == 'NA':   # if source is not given
		for node in G.nodes():
			if G.out_degree(node) == 0:
				source = node
		print('root', source)

	# performing bfs
	G2 = G.to_undirected()   # convert graph to undirected
	T = nx.bfs_tree(G2, source=source)
	# print(list(T.nodes()))
	# print(list(T.edges()))

	# assign qs number
	assigndict = {list(T.nodes())[0]: 1}   # assign the first node 
	# print(assigndict)
	for node in list(T.nodes()):
		# print('node', node)
		# get all neighbors, including unassigned and assigned neighbors
		nbs = list(G2.neighbors(node))
		assigned_nbs = [nb for nb in nbs if nb in assigndict.keys()]
		unassigned_nbs = Diff(nbs, assigned_nbs)
		# print('assigned nbs', assigned_nbs)
		# print('unassigned nbs', unassigned_nbs)
		# get available colors (colors unsed by itself or neighboring cells)
		used_color = [assigndict[n] for n in assigned_nbs+[node] if n in assigndict.keys()]
		# print('used color', used_color)
		avail_color = Diff(list(colorcode), used_color)
		# print(avail_color)
		# assign available color to unassigned neighbor (random assignment)
		for nb in unassigned_nbs:
			assigndict[nb] = random.choice(avail_color) 
			avail_color.remove(assigndict[nb])   # remove this color from available colors
		# print(assigndict)
	return T, assigndict


def place_assignment_on_hexmap (G, assigndict, T, hexmap):
	# generate a color dictionary
	color = sns.color_palette('hls', n_colors=4)
	colordict = {}
	colorcode = np.arange(1, 5)
	for k in colorcode:
		colordict[k] = color[k-1]

	placement_tmp = copy.deepcopy (G)
	G2 = G.to_undirected()   # convert graph to undirected
	# randomly choose a node from hexagon map that has a colorcode of 1
	colorcode_1 = [node for node in hexmap.nodes() if hexmap.nodes[node]['colorcode'] == 1]
	boundary = max([hexmap.nodes[node]['coor'][0] for node in hexmap.nodes()]) 
	center_locs = [node for node in colorcode_1 if (hexmap.nodes[node]['coor'][0]<=0.75*boundary and hexmap.nodes[node]['coor'][0]>=0.25*boundary) and (hexmap.nodes[node]['coor'][1]<=0.75*boundary and hexmap.nodes[node]['coor'][1]>=0.25*boundary)]
	# print(center_locs)
	fp = random.choice(center_locs)   # place the root cell in the randomly selected hexmap node
	# print('placing root at', fp)
	placement_tmp.nodes[list(T.nodes())[0]]['placement'] = fp
	placed_cells = [list(T.nodes())[0]]
	for cell in list(T.nodes()):
		nbs = list(G2.neighbors(cell))
		for nb in nbs:
			if nb not in placed_cells:
				# print(nb,'not placed')
				avail_locs = list(hexmap.neighbors(placement_tmp.nodes[cell]['placement']))
				# print('available locs', avail_locs)
				placement_tmp.nodes[nb]['placement'] = [loc for loc in avail_locs if hexmap.nodes[loc]['colorcode'] == assigndict[nb]][0]
				# print('placing node at loc', placement_tmp.nodes[nb]['placement'])
				placed_cells.append(nb)
	return placement_tmp

def find_allowed_placement_on_hexmap (G, hexmap, source, steps):
	""" find a color assignment that has no two nodes occupying the same location on hexagon map """
	allowed = False
	trial = 0
	while allowed == False:
		T, assigndict = assign_colors (G, source)
		placement_tmp = place_assignment_on_hexmap (G, assigndict, T, hexmap)
		locs = [placement_tmp.nodes[cell]['placement'] for cell in placement_tmp.nodes()]
		if len(locs) == len(set(locs)): 
			print('found placement solution')
			allowed = True 
		trial += 1
		if trial > steps: 
			print('no placement solution found, adding relayers')
			print('T', T)
			for dist in range(T):
				print(dist, nx.descendants_at_distance(G, source, dist))
			break 
	return placement_tmp, allowed


def place_assignment_on_hexmap_wbuffer (G, hexmap, source, steps):
	# generate a color dictionary
	color = sns.color_palette('hls', n_colors=4)
	colordict = {}
	colorcode = np.arange(1, 5)
	for k in colorcode:
		colordict[k] = color[k-1]

	# start from a node that has no outgoing edge
	if source == 'NA':   # if source is not given
		for node in G.nodes():
			if G.out_degree(node) == 0:
				source = node
		print('root', source)

	# performing bfs
	G2 = G.to_undirected()   # convert graph to undirected
	T = nx.bfs_tree(G2, source=source)
	print('bfs nodes', list(T.nodes()))
	print('bfs edges', list(T.edges()))

	allowed = False
	trial = 0

	while allowed == False:

		# place cells on grid 
		placement_tmp = copy.deepcopy (G2)
		buffer_num = 0

		# randomly choose a node from hexagon map that has a colorcode of 1
		colorcode_1 = [node for node in hexmap.nodes() if hexmap.nodes[node]['colorcode'] == 1]
		boundary = max([hexmap.nodes[node]['coor'][0] for node in hexmap.nodes()]) 
		center_locs = [node for node in colorcode_1 if (hexmap.nodes[node]['coor'][0]<=0.6*boundary and hexmap.nodes[node]['coor'][0]>=0.4*boundary) and (hexmap.nodes[node]['coor'][1]<=0.6*boundary and hexmap.nodes[node]['coor'][1]>=0.4*boundary)]
		# print(center_locs)
		fp = random.choice(center_locs)   # place the root cell in the randomly selected hexmap node
		print('placing root at', fp)
		placement_tmp.nodes[list(T.nodes())[0]]['placement'] = fp
		placed_cells = [list(T.nodes())[0]]
		occupied_hexlocs = [fp]

		try:
			for cell in list(T.nodes()):
				print('placing neighbor cells of ', cell)
				nbs = list(set(list(G2.neighbors(cell))) - set(placed_cells))  # not placed neighbors
				print('neighbors', nbs)
				avail_locs = list(set(list(hexmap.neighbors(placement_tmp.nodes[cell]['placement']))) - set(occupied_hexlocs))
				print('available locs', avail_locs)
				if len(avail_locs) >= len(nbs):   # if there is enough space to place the cell
					for nb in nbs:
						print(nb,'not placed')
						placement = random.choice(avail_locs)
						placement_tmp.nodes[nb]['placement'] = placement
						print('placing node at loc', placement)
						avail_locs.remove(placement)
						occupied_hexlocs.append(placement)
						placed_cells.append(nb)
				else:   # if not enough space to place the cell, add buffer cell
					if len(avail_locs) != 0:
						print("adding buffer cell")
						# add buffer cell
						buffer_cell = 'B'+str(buffer_num+1)
						buffer_num += 1
						placement_tmp.add_node(buffer_cell)
						placement_tmp.add_edge(cell, buffer_cell)
						for nb in nbs:
							placement_tmp.add_edge(buffer_cell, nb)
							placement_tmp.remove_edge(cell, nb)
						# place buffer cell
						placement = random.choice(avail_locs)
						print('placing buffer cell', buffer_cell, 'at', placement)
						placement_tmp.nodes[buffer_cell]['placement'] = placement
						occupied_hexlocs.append(placement)
						buf_avail_locs = list(set(list(hexmap.neighbors(placement))) - set(occupied_hexlocs))
						for nb in nbs: 
							print(nb, 'not placed')
							placement = random.choice(buf_avail_locs)
							placement_tmp.nodes[nb]['placement'] = placement
							print('placing node at loc', placement)
							buf_avail_locs.remove(placement)
							occupied_hexlocs.append(placement)
							placed_cells.append(nb)
					else: 
						print('no solution in trial', trial)
		except KeyError: 
			print('no solution in trial', trial)

		trial += 1
		# check whether there's overlap in occupied positions on the hexmap
		locs = [placement_tmp.nodes[cell]['placement'] for cell in placement_tmp.nodes()]
		if len(locs) == len(set(locs)) and len(locs) == len(list(placement_tmp.nodes())): 
			print('found placement solution')
			allowed = True 
			break
		if trial > steps: 
			break 

	return placement_tmp, allowed


##########################################################
# Visualize placement result
##########################################################

def visualize_hexmap_placement_graphviz (hexmap, G_placement, outfile):
	""" visualize graph placement on hexagonal grid """

	# create color dictionary
	color = sns.color_palette('hls', n_colors=4)
	colordict = {}
	colorcode = np.arange(1, 5)
	for k in colorcode:
		colordict[k] = color[k-1]

	# draw grid nodes
	X, Y = [], []
	for node in hexmap.nodes():
		coor = hexmap.nodes[node]['coor']
		plt.plot([coor[0]], [coor[1]], markersize=1, marker='o', color='lightgrey')

	# draw grid edges
	for edge in hexmap.edges():
		n0, n1 = edge[0], edge[1]
		c0, c1 = hexmap.nodes[n0]['coor'], hexmap.nodes[n1]['coor']
		xlist = [c0[0], c1[0]]
		ylist = [c0[1], c1[1]]
		plt.plot(xlist, ylist, 'lightgrey', zorder=0)

	# draw circuit gates
	for node in G_placement.nodes():
		placement = G_placement.nodes[node]['placement'] 
		coor = hexmap.nodes[placement]['coor']
		plt.plot([coor[0]], [coor[1]], markersize=2, marker='o', color=colordict[hexmap.nodes[placement]['colorcode']])
		plt.text(coor[0], coor[1], node)

	# draw circuit edges
	for edge in G_placement.edges():
		n0, n1 = edge[0], edge[1]
		c0, c1 = hexmap.nodes[G_placement.nodes[n0]['placement']]['coor'], hexmap.nodes[G_placement.nodes[n1]['placement']]['coor']
		xlist = [c0[0], c1[0]]
		ylist = [c0[1], c1[1]]
		plt.plot(xlist, ylist, 'k', zorder=0)

	plt.box(on=None)
	plt.savefig(outfile+'.pdf', dpi=200)
	plt.show()	



if __name__ == '__main__':

	PATH = '/Users/jgzhang/Work/Densmore_lab/genetic-circuit-partitioning/2021.4/runs/benchmark/examples-placing/tree3'

	hexmap = initialize_hexmap (20)
	assign_hexmap_color (hexmap, 20)
	# visualize_hexmap (hexmap, './examples/tree3/grid_20')

	# path = './examples/374'
	# G, innodes, outnodes = load_graph (path+'/DAG.edgelist')
	# G_p = return_prim_graph (G, innodes, outnodes)
	# visualize_graph_graphviz (G_p, './examples/374_dag')

	# G = generate_graph (20, PATH+'/examples/tree3')
	G = load_undirected_graph (PATH+'/tree3.edgelist')
	# visualize_graph_graphviz (G, './examples/tree3_DAG')
	# T, assigndict = assign_colors (G, 1)
	# G_placement = place_assignment_on_hexmap (G, assigndict, T, hexmap)
	# visualize_hexmap_placement_graphviz (hexmap, G_placement, './examples/random_tree')
	# G_placement, allowed = find_allowed_placement_on_hexmap (G, hexmap, 1, 1000)
	G_placement, allowed = place_assignment_on_hexmap_wbuffer (G, hexmap, 1, 10)
	print('Found solution', allowed)
	if allowed: 
		visualize_hexmap_placement_graphviz (hexmap, G_placement, PATH+'/tree3')

