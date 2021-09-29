import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import genetic_partition_test as gp

def load_graph (edgelist):
    G = nx.read_edgelist (edgelist, nodetype = str, create_using=nx.DiGraph())
    return G

def community_layout(g, partition):
    """
    Compute the layout for a modular graph.


    Arguments:
    ----------
    g -- networkx.Graph or networkx.DiGraph instance
        graph to plot

    partition -- dict mapping int node -> int community
        graph partitions


    Returns:
    --------
    pos -- dict mapping int node -> (float x, float y)
        node positions

    """

    pos_communities = _position_communities(g, partition, scale=3.)

    pos_nodes = _position_nodes(g, partition, scale=1.)

    # combine positions
    pos = dict()
    for node in g.nodes():
        print ('community position', pos_communities[node])
        print ('node position', pos_nodes[node])
        pos[node] = pos_communities[node] + pos_nodes[node]
        print('new node position', pos[node])

    return pos

def _position_communities(g, partition, **kwargs):

    # create a weighted graph, in which each node corresponds to a community,
    # and each edge weight to the number of edges between communities
    between_community_edges = _find_between_community_edges(g, partition)

    communities = set(partition.values())
    hypergraph = nx.DiGraph()
    hypergraph.add_nodes_from(communities)
    for (ci, cj), edges in between_community_edges.items():
        hypergraph.add_edge(ci, cj, weight=len(edges))

    # find layout for communities
    pos_communities = nx.spring_layout(hypergraph, **kwargs)

    # set node positions to position of community
    pos = dict()
    for node, community in partition.items():
        pos[node] = pos_communities[community]

    return pos

def _find_between_community_edges(g, partition):

    edges = dict()

    for (ni, nj) in g.edges():
        ci = partition[ni]
        cj = partition[nj]

        if ci != cj:
            try:
                edges[(ci, cj)] += [(ni, nj)]
            except KeyError:
                edges[(ci, cj)] = [(ni, nj)]

    return edges

def _position_nodes(g, partition, **kwargs):
    """
    Positions nodes within communities.
    """

    communities = dict()
    for node, community in partition.items():
        try:
            communities[community] += [node]
        except KeyError:
            communities[community] = [node]

    pos = dict()
    for ci, nodes in communities.items():
        subgraph = g.subgraph(nodes)
        pos_subgraph = nx.spring_layout(subgraph, **kwargs)
        pos.update(pos_subgraph)

    return pos

def plot_partition ():

    outdir = "/Users/jgzhang/Work/Densmore_lab/Partition/code_version/v2/genetic-circuit-partitioning/2021.4/"
    bm_path = outdir + 'runs/benchmark/electronic-circuits/RCA4/'
    sol_path = outdir + 'runs/results/electronic-circuits/RCA4/nparts/14/optimized_lc/'
    part_sol = sol_path + 'part_solns.txt'

    edgelist = bm_path + '/DAG.edgelist'
    G = load_graph (edgelist)
    in_nodes, out_nodes, nonprimitives  = gp.get_nonprimitive_nodes (G)
    G_primitive = gp.get_G_primitive (G, nonprimitives)

    solDict = gp.load_opt_part_sol (part_sol)
    iteration = 6

    part = solDict[iteration]['part']
    node_to_partDict = {}
    for n in G_primitive.nodes():
        node_to_partDict[n] = gp.get_part(part, n)

    pos = community_layout(G_primitive, node_to_partDict)

    nx.draw(G_primitive, pos, node_color=list(node_to_partDict.values()))
    plt.show()
    return


if __name__ == '__main__':
    plot_partition() 