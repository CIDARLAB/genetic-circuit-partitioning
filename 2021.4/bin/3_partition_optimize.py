#!/usr/bin/env python

# 	Copyright (C) 2021 by
# 	Jing Zhang <jgzhang@bu.edu>, Densmore Lab, Boston University
# 	All rights reserved.
# 	OSI Non-Profit Open Software License ("Non-Profit OSL") 3.0 license.

# Supporting modules
import argparse
import genetic_partition as gp
import os


def main():

    # Parse the command line inputs
    parser = argparse.ArgumentParser(description="perform graph partition using metis")
    parser.add_argument(
        "-settings",
        dest="settings",
        required=True,
        help="settings.txt",
        metavar="string",
    )
    parser.add_argument(
        "-samples", dest="samples", required=True, help="1,2", metavar="string"
    )
    args = parser.parse_args()

    # Run the command
    samples = args.samples.split(",")
    settings = gp.load_settings(args.settings)
    for s in samples:
        print("Processing sample", s)
        # obtain user-defined params
        S_bounds = settings[s]["S_bounds"].split(",")
        target_n = settings[s]["target_n"].split(",")
        primitive_only = settings[s]["primitive_only"]
        maxNodes = int(settings[s]["max_nodes_to_shift"])
        high_constraint = settings[s]["high_constraint"].split(",")
        low_constraint = settings[s]["low_constraint"].split(",")
        trajectories = int(settings[s]["trajectories"])
        out_path = settings[s]["output_path"]

        # load graph
        G = gp.load_graph_undirected(settings, s)
        DAG = gp.load_graph(settings, s)

        outdir = out_path + "/nparts/"
        order = gp.rank_connectivity(DAG, primitive_only, outdir)
        print("median degree of connectivity in subgraphs", order)

        # optimize graph partition results to satisfy constraints
        if target_n == [""]:
            for npart in os.listdir(outdir):
                if npart.isdigit():
                    print("target npart", npart)
                    part_sol = outdir + npart + "/part_solns.txt"
                    cut, partDict = gp.load_metis_part_sol(part_sol)
                    if len(list(G.nodes())) >= 25:
                        print("optimizing by subnetworks")
                        gp.optimize_signal_subnetwork(
                            DAG,
                            primitive_only,
                            S_bounds,
                            cut,
                            partDict,
                            maxNodes,
                            low_constraint,
                            trajectories,
                            outdir + npart + "/optimized_lc/",
                        )
                        gp.optimize_signal_subnetwork(
                            DAG,
                            primitive_only,
                            S_bounds,
                            cut,
                            partDict,
                            maxNodes,
                            high_constraint,
                            trajectories,
                            outdir + npart + "/optimized_hc/",
                        )
                    else:
                        print("brute-force optimizing")
                        gp.optimize_signal_bruteforce(
                            DAG,
                            primitive_only,
                            S_bounds,
                            cut,
                            partDict,
                            maxNodes,
                            low_constraint,
                            outdir + npart + "/optimized_lc/",
                        )
                        gp.optimize_signal_bruteforce(
                            DAG,
                            primitive_only,
                            S_bounds,
                            cut,
                            partDict,
                            maxNodes,
                            low_constraint,
                            outdir + npart + "/optimized_hc/",
                        )
            gp.determine_best_solution(
                DAG, primitive_only, high_constraint, low_constraint, outdir
            )

        else:
            for n in target_n:
                print("target npart", n)
                outdir = out_path + "/nparts/"
                part_sol = outdir + n + "/part_solns.txt"
                cut, partDict = gp.load_metis_part_sol(part_sol)
                if len(list(G.nodes())) >= 25:
                    print("optimizing by subnetworks")
                    gp.optimize_signal_subnetwork(
                        DAG,
                        primitive_only,
                        S_bounds,
                        cut,
                        partDict,
                        maxNodes,
                        low_constraint,
                        trajectories,
                        outdir + n + "/optimized_lc/",
                    )
                    gp.optimize_signal_subnetwork(
                        DAG,
                        primitive_only,
                        S_bounds,
                        cut,
                        partDict,
                        maxNodes,
                        high_constraint,
                        trajectories,
                        outdir + n + "/optimized_hc/",
                    )
                else:
                    print("brute-force optimizing")
                    gp.optimize_signal_bruteforce(
                        DAG,
                        primitive_only,
                        S_bounds,
                        cut,
                        partDict,
                        maxNodes,
                        low_constraint,
                        outdir + n + "/optimized_lc/",
                    )
                    gp.optimize_signal_bruteforce(
                        DAG,
                        primitive_only,
                        S_bounds,
                        cut,
                        partDict,
                        maxNodes,
                        low_constraint,
                        outdir + n + "/optimized_hc/",
                    )
            gp.determine_best_solution(
                DAG, primitive_only, high_constraint, low_constraint, outdir
            )


if __name__ == "__main__":
    main()
