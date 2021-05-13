import random
import matplotlib.pyplot as plt
from matplotlib import cm
import networkx as nx
import metis
from collections import Counter
import numpy as np
import time
import numpy.linalg as la
import scipy.cluster.vq as vq
import itertools
import operator
import math
import copy
from mpmath import *
from itertools import chain
from itertools import product
from itertools import starmap
from functools import partial
import os
from networkx.drawing.nx_agraph import graphviz_layout


def load_gate_info(inputfile):
    lines = [open(inputfile, "r").read().strip("\n")][0].split("\n")
    gateLib = {}
    header = lines[0].split("\t")
    for line in lines[1:]:
        tokens = line.split("\t")
        gate = tokens[0]
        gateLib[gate] = {}
        for idx, token in enumerate(tokens, 1):
            gateLib[gate][header[idx - 1]] = token
    return gateLib


def load_sensor_info(inputfile):
    lines = [open(inputfile, "r").read().strip("\n")][0].split("\n")
    sensorLib = {}
    for line in lines[1:]:
        tokens = line.split("\t")
        sensor = tokens[0]
        REU_off = float(tokens[1])
        REU_on = float(tokens[2])
        sensorLib[sensor] = {"on": REU_on, "off": REU_off}
    return sensorLib


def sigmoid(ymin, ymax, Kd, n, x):
    return ymin + (ymax - ymin) / (1 + math.pow((x / Kd), n))


def sigmoid_curve(xlist, ymin, ymax, Kd, n):
    ylist = [sigmoid(ymin, ymax, Kd, n, x) for x in xlist]
    return ylist


def load_logic(inputfile):
    lines = [open(inputfile, "r").read().strip("\n")][0].split("\n")
    logic = {"PTac": [], "PTet": [], "PBAD": []}
    for line in lines[1:]:
        tokens = line.split("\t")[1].split(",")
        logic["PTac"].append(int(tokens[0]))
        logic["PTet"].append(int(tokens[1]))
        logic["PBAD"].append(int(tokens[2]))
    return logic


def load_logic_4input(inputfile):
    lines = [open(inputfile, "r").read().strip("\n")][0].split("\n")
    logic = {"PTac": [], "PTet": [], "PBAD": [], "PCin": []}
    for line in lines[1:]:
        tokens = line.split("\t")[1].split(",")
        logic["PTac"].append(int(tokens[0]))
        logic["PTet"].append(int(tokens[1]))
        logic["PBAD"].append(int(tokens[2]))
        logic["PCin"].append(int(tokens[3]))
    return logic


def calc_circuit_score(assignedLib, gateLib, sensorLib):
    logic = load_logic_4input("./gate assignment/statelogic-4input.txt")
    circuit_valid = True

    # initialize a score dictionary
    scoreDict = {}
    for node in assignedLib:
        if assignedLib[node]["type"] != "input":
            scoreDict[node] = {"logic": [], "output REU": []}

    for i in range(16):
        visited = 0
        assignedLib_i = copy.deepcopy(assignedLib)  # set up a temporary library

        # first calculate the input REU
        for node in assignedLib_i:
            if assignedLib_i[node]["type"] == "input":
                if logic[assignedLib_i[node]["sensor"]][i] == 0:
                    assignedLib_i[node]["output REU"] = assignedLib_i[node]["REU OFF"]
                    assignedLib_i[node]["logic"] = 0
                else:
                    assignedLib_i[node]["output REU"] = assignedLib_i[node]["REU ON"]
                    assignedLib_i[node]["logic"] = 1
                # print(i, node, assignedLib_i[node]['sensor'], assignedLib_i[node])
                visited += 1

        # calculate the REU of primitive and output node
        for node in assignedLib_i:
            if len(assignedLib_i[node]["in"]) == 1:
                assignedLib_i[node]["visited"] = [0]
                assignedLib_i[node]["logic"] = [-1]
            elif len(assignedLib_i[node]["in"]) == 2:
                assignedLib_i[node]["visited"] = [0, 0]
                assignedLib_i[node]["logic"] = [-1, -1]
        r = 1
        while visited != len(assignedLib_i.keys()):
            # print('round ###########################################', r)
            for node in assignedLib_i:
                if assignedLib_i[node]["output REU"] == 0:
                    # print('node', node)
                    # get input REU
                    # print('incoming nodes', assignedLib_i[node]['in'])
                    in_nodes = assignedLib_i[node]["in"]
                    for idx, in_node in enumerate(in_nodes):
                        # print('input node', in_node, assignedLib_i[in_node]['output REU'], assignedLib_i[in_node]['logic'])
                        # if the in node has a calculated output REU
                        if (
                            assignedLib_i[in_node]["output REU"] != 0
                            and assignedLib_i[node]["visited"][idx] != 1
                        ):
                            # print(in_node, assignedLib_i[in_node]['output REU'])
                            assignedLib_i[node]["input REU"].append(
                                assignedLib_i[in_node]["output REU"]
                            )
                            assignedLib_i[node]["visited"][idx] = 1
                            assignedLib_i[node]["logic"][idx] = assignedLib_i[in_node][
                                "logic"
                            ]

                    # print('inputREU', assignedLib_i[node]['input REU'])
                    # print(assignedLib_i[node]['visited'])
                    # output REU
                    if 0 not in assignedLib_i[node]["visited"]:
                        if assignedLib[node]["type"] != "output":
                            params = assignedLib_i[node]["params"]
                            x = sum(assignedLib_i[node]["input REU"])
                            # print('x', x)
                            # print('params', params)
                            assignedLib_i[node]["output REU"] = sigmoid(
                                params[0], params[1], params[2], params[3], x
                            )
                            if 1 in assignedLib_i[node]["logic"]:
                                assignedLib_i[node]["logic"] = 0
                            else:
                                assignedLib_i[node]["logic"] = 1
                            # print('output REU', assignedLib_i[node]['output REU'])
                            # print('logic of gate', node, assignedLib_i[node]['logic'])
                            # print('number of gates that have output REU', visited)
                        else:
                            assignedLib_i[node]["output REU"] = sum(
                                assignedLib_i[node]["input REU"]
                            )
                            if 1 in assignedLib_i[node]["logic"]:
                                assignedLib_i[node]["logic"] = 1
                            else:
                                assignedLib_i[node]["logic"] = 0
                            # print('done')
                            # update score dictionary

                        scoreDict[node]["logic"].append(assignedLib_i[node]["logic"])
                        scoreDict[node]["output REU"].append(
                            assignedLib_i[node]["output REU"]
                        )
                        visited += 1
            r += 1

    # calculate score of this permutation
    Smin = 1e3
    for node in assignedLib:
        if assignedLib[node]["type"] == "output":
            # print(node)
            # print(scoreDict[node]['output REU'])
            # print(scoreDict[node]['logic'])
            try:
                maxOFF = max(
                    [
                        scoreDict[node]["output REU"][s]
                        for s in range(8)
                        if scoreDict[node]["logic"][s] == 0
                    ]
                )
                minON = min(
                    [
                        scoreDict[node]["output REU"][s]
                        for s in range(8)
                        if scoreDict[node]["logic"][s] == 1
                    ]
                )
                # print('min on', minON, 'max off', maxOFF)
                S = minON / maxOFF
                if S < Smin:
                    Smin = S
            except ValueError:
                scoreDict[node]["logic"]
                circuit_valid = False
    # return the lowest S
    return Smin, circuit_valid


def record_library(assignedLib, outfile):
    """ record the library with gate assignments that generate the highest circuit score """
    f_out = open(outfile, "w")
    f_out.write("node\ttype\tsensor/gate\tparams\tpartition\n")
    for node in assignedLib:
        if assignedLib[node]["type"] == "input":
            f_out.write(
                "\t".join(
                    [
                        node,
                        assignedLib[node]["type"],
                        assignedLib[node]["sensor"],
                        str(
                            [assignedLib[node]["REU ON"], assignedLib[node]["REU OFF"]]
                        ),
                        "na",
                    ]
                )
                + "\n"
            )
        elif assignedLib[node]["type"] == "primitive":
            f_out.write(
                "\t".join(
                    [
                        node,
                        assignedLib[node]["type"],
                        assignedLib[node]["gate"],
                        str(assignedLib[node]["params"]),
                        str(assignedLib[node]["part"]),
                    ]
                )
                + "\n"
            )
        else:
            f_out.write(
                "\t".join([node, assignedLib[node]["type"], "na", "na", "na"]) + "\n"
            )


def assign_gates(G, edge_direct, partition, nonprimitives, params, n, k):
    # assign biological gates to partitioned graphs
    gateLib = load_gate_info("./gate assignment/gatefunction.txt")
    sensorLib = load_sensor_info("./gate assignment/sensor.txt")
    # print(gateLib)

    assignedLib = {}
    # add 'input node' and 'output node' of each gate
    for v in G.nodes():
        assignedLib[v] = {}
        assignedLib[v]["in"], assignedLib[v]["out"] = [], []
        assignedLib[v]["input REU"] = []
        assignedLib[v]["output REU"] = 0
        print(v, assignedLib[v])
    for e in edge_direct:
        print(e)
        assignedLib[e[0]]["out"].append(e[1])
        assignedLib[e[1]]["in"].append(e[0])

    # assign input nodes
    # for ADDER
    # inodes = ['x', 'y', 'cin']
    # onodes = ['A', 'cout']
    # for synthesize graphs
    inodes = ["a0", "a1", "b0", "b1"]
    onodes = ["o0", "o1", "o2", "o3"]
    sensors = ["PTac", "PTet", "PBAD", "PCin"]
    for idx, v in enumerate(inodes):
        assignedLib[v]["REU ON"] = sensorLib[sensors[idx]]["on"]
        assignedLib[v]["REU OFF"] = sensorLib[sensors[idx]]["off"]
        assignedLib[v]["sensor"] = sensors[idx]
        assignedLib[v]["type"] = "input"
        print(assignedLib[v])

    for v in onodes:
        assignedLib[v]["type"] = "output"

    # copy another G, remove input and output nodes
    G_primitive = copy.deepcopy(G)
    for node in nonprimitives:
        G_primitive.remove_node(node)

    ## initialize the gate assignment
    for i in range(0, max(partition[1]) + 1):
        # print(partition[1])
        nodeIdx = [a for a, b in enumerate(partition[1]) if b == i]
        gnodes = [list(G_primitive.nodes())[n] for n in nodeIdx]
        # print(gnodes, 'gate nodes within this partition')
        # update partition result
        for v in gnodes:
            assignedLib[v]["part"] = i  # partition of this gate
        # assign repressor gate
        chosengates = random.sample(range(1, 13), len(gnodes))
        # print(chosengates)
        for idx, v in enumerate(gnodes):
            vg = chosengates[idx]
            vg_rbs = [
                key for key in gateLib if key.split("-")[0] == str(vg)
            ]  # different rbs choices for this repressor gate
            vg = random.choice(vg_rbs)
            assignedLib[v]["gate"] = vg
            assignedLib[v]["params"] = [
                float(gateLib[vg]["ymin"]),
                float(gateLib[vg]["ymax"]),
                float(gateLib[vg]["K"]),
                float(gateLib[vg]["n"]),
            ]
            assignedLib[v]["type"] = "primitive"

        # calculate circuit score of initial assignments
    S0, circuit_valid = calc_circuit_score(assignedLib, gateLib, sensorLib)

    if circuit_valid:
        # simulated annealing
        Tmax = params[0]  # starting temperature
        C = params[1]

        for i in range(5):  # for 100 trajectories
            SList = []
            Smax = S0
            bestLib = copy.deepcopy(assignedLib)

            for t in range(1000000):  # perform simulated annealing
                # print('trial', t)
                # clone a assignedLib
                assignedLib_tmp = copy.deepcopy(assignedLib)
                # print('new tmp assignedLib', assignedLib_tmp)
                # randomly select a node
                g_swap = random.choice(list(G_primitive.nodes()))
                g_part = assignedLib_tmp[g_swap]["part"]
                # print(g_swap, g_part, assignedLib_tmp[g_swap]['gate'], 'gate to be swapped')
                # get all available gates (a gate already in this partition, a different gate with a never used repressor from the library)
                # print('gates within the gate library', list(gateLib.keys()))
                used_gates = [
                    assignedLib_tmp[g]["gate"]
                    for g in assignedLib_tmp
                    if assignedLib_tmp[g]["type"] == "primitive"
                    and assignedLib_tmp[g]["part"] == g_part
                    and g != g_swap
                ]
                # print('other used gate within this partition', used_gates)
                used_repr = [r.split("-")[0] for r in used_gates]
                # print('other used repressor within this partition', used_repr)
                # remove repressors that have already been used by other nodes in this partition
                availgates = list(
                    set(gateLib.keys())
                    - set([g for g in gateLib if g.split("-")[0] in used_repr])
                )
                availgates.remove(assignedLib_tmp[g_swap]["gate"])
                # add other gates that have been used in this partition
                availgates = availgates + used_gates
                # print('available gates', availgates)

                # swap two gates g_swap and g_swap_f
                g_swap_t = random.choice(availgates)
                # print('gate swapped to', g_swap_t)
                if (
                    g_swap_t in used_gates
                ):  # if swapping two gates within the same partition
                    # update the gate and transfer function params of the swapped gate
                    g_swap_f = [
                        g
                        for g in assignedLib_tmp
                        if assignedLib_tmp[g]["type"] == "primitive"
                        and assignedLib_tmp[g]["part"] == g_part
                        and assignedLib_tmp[g]["gate"] == g_swap_t
                    ][0]
                    assignedLib_tmp[g_swap_f]["gate"] = assignedLib_tmp[g_swap]["gate"]
                    assignedLib_tmp[g_swap_f]["params"] = [
                        float(gateLib[assignedLib_tmp[g_swap]["gate"]]["ymin"]),
                        float(gateLib[assignedLib_tmp[g_swap]["gate"]]["ymax"]),
                        float(gateLib[assignedLib_tmp[g_swap]["gate"]]["K"]),
                        float(gateLib[assignedLib_tmp[g_swap]["gate"]]["n"]),
                    ]
                    # update the transfer function params of the chosen gate
                    assignedLib_tmp[g_swap]["gate"] = g_swap_t
                    assignedLib_tmp[g_swap]["params"] = [
                        float(gateLib[g_swap_t]["ymin"]),
                        float(gateLib[g_swap_t]["ymax"]),
                        float(gateLib[g_swap_t]["K"]),
                        float(gateLib[g_swap_t]["n"]),
                    ]
                    # print(assignedLib[g_swap])
                    # print(assignedLib_tmp[g_swap])
                    # print(assignedLib[g_swap_f])
                    # print(assignedLib_tmp[g_swap_f])
                else:  # if swapping with a gate in the gate library
                    assignedLib_tmp[g_swap]["gate"] = g_swap_t
                    # update the transfer function params of the chosen gate
                    assignedLib_tmp[g_swap]["params"] = [
                        float(gateLib[g_swap_t]["ymin"]),
                        float(gateLib[g_swap_t]["ymax"]),
                        float(gateLib[g_swap_t]["K"]),
                        float(gateLib[g_swap_t]["n"]),
                    ]
                    # print(assignedLib[g_swap])
                    # print(assignedLib_tmp[g_swap])

                # calculate the S score, store it
                S, circuit_valid = calc_circuit_score(
                    assignedLib_tmp, gateLib, sensorLib
                )
                # print('score', S, 'original score', S0)

                # choose to accept or reject the choice

                Ti = Tmax * (math.exp(-C * t))
                if Ti > 0:
                    try:
                        P = math.exp((S - S0) / Ti)
                    except OverflowError:
                        P = 2
                    # print('temperature', Ti)
                    # print('Probability', P)
                else:
                    P = math.exp((S - S0) / 0.01)

                # append the highest score that occurs to this point
                if S > Smax:
                    Smax = S
                    bestLib = copy.deepcopy(assignedLib_tmp)
                SList.append(Smax)
                # print('highest score that occurs to round', t, SList)

                # if P > 1, accept change
                if P > 1:
                    assignedLib = assignedLib_tmp
                    # print('P>1, accepted')
                # if P < 1, accept change based on probability
                else:
                    if random.random() < P:
                        assignedLib = assignedLib_tmp
                        # print('P<1, accepted')
                    # else:
                    # print('P<1, rejected')

                # print('assigned library after round', t, assignedLib)

            # record the max S
            # path = './Random graphs/n'+str(n)+'/DAG'+str(n)+'.'+str(k)
            # if not os.path.exists(path):
            # 	os.mkdir(path)
            path = "/Users/jgzhang/Work/Densmore lab/Partition/COI/4input-305"
            try:
                outfile = open(path + "/S_trajectory.txt", "a")
            except FileNotFoundError:
                outfile = open(path + "/S_trajectory.txt", "w")

            outfile.write(str(i) + "\t")
            outfile.write(",".join([str(S) for S in SList]) + "\n")
            # record the best library
            record_library(
                bestLib,
                path + "/optimal gate assignments_trajectory_" + str(i) + ".txt",
            )

    return G


def plot_circuit_score(inputfile, outfile):
    """ plot circuit score trajectories in simulated annealing """
    lines = [open(inputfile, "r").read().strip("\n")][0].split("\n")
    fig = plt.figure(figsize=(5, 3))
    ax = fig.add_subplot(111)
    for line in lines:
        SList = [float(s) for s in line.split("\t")[1].split(",")]
        print(line[0], SList[-1])
        ax.plot(range(len(SList)), SList, linewidth=0.2)
    # plt.ylim([0, 350])
    plt.savefig(outfile, dpi=200)
    plt.show()


def load_gate_assignment(S_file, G):
    """ load partitioned circuit with biological gates assigned """
    # lines = [open(S_file), 'r').read().strip("\n")][0].split('\n')
    # get the best run with the highest circuit score
    bestrun, bestS = "NA", 0
    for line in lines:
        run, SList = line.split("\t")[0], [
            float(s) for s in line.split("\t")[1].split(",")
        ]
        if SList[-1] > bestS:
            bestrun = run
    # load the gate assignment data from the best run
    lines = [
        open(opt_gate_assignment_filename(settings, sample, bestrun), "r")
        .read()
        .strip("\n")
    ][0].split("\n")

    # update G with gate assignment information
    for line in lines[1:]:
        tokens = line.split("\t")
        node = tokens[0]
        G.nodes[node]["type"] = tokens[1]
        G.nodes[node]["gate"] = tokens[2]
        G.nodes[node]["params"] = tokens[3]
        G.nodes[node]["part"] = tokens[4]
    return G


def assign_gates_coi(G, partition, nonprimitives):
    toxicity = {
        "AmeR": 1.06,
        "AmtR": 1.29,
        "BetI": 0.94,
        "BM3R1": 1.01,
        "HlyIIR": 1.16,
        "IcaRA": 1.93,
        "LitR": 1.48,
        "LmrA": 1,
        "PhlF": 0.97,
        "PsrA": 1,
        "QacR": 1.95,
        "SrpR": 1.07,
    }
    for i in range(0, max(partition[1]) + 1):
        avail_gates = list(toxicity.keys())
        nodeIdx = [a for a, b in enumerate(partition[1]) if b == i]
        nodes = [list(G.nodes())[n] for n in nodeIdx]
        print(nodes)
        for v in nodes:
            if v not in nonprimitives:
                # if the node is not input or output
                print(v)
                G.nodes[v]["gate"] = random.choice(avail_gates)
                G.nodes[v]["weight"] = toxicity[G.nodes[v]["gate"]]
                avail_gates.remove(G.nodes[v]["gate"])
            else:
                # add a lot of weight to input or output
                G.nodes[v]["weight"] = 99
    return G
