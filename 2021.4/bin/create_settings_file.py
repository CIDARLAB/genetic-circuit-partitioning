import os


PATH = "/Users/jgzhang/Work/Densmore_lab/genetic-circuit-partitioning/2021.4/runs"
output = PATH + "/5-input-settings.txt"
benchmarks = PATH + "/benchmark/5-input-boolean-circuits/"
graphs = ""

f_out = open(output, "w")

for graph in os.listdir(benchmarks):
    if graph.isdigit():
        # print(graph)
        graphs = graphs + "," + graph
        # print(graphs)
        graph_path = "./benchmark/5-input-boolean-circuits/" + graph
        lib_path = "./lib"
        primitive_only = "True"
        S_bounds = "3,7"
        target_n = ""
        maxNodes = "5"
        hc = "2,2"
        lc = "4"
        t = "10"
        output_path = "./results/5-input-boolean-circuits/" + graph
        f_out.write(
            "\t".join(
                [
                    graph,
                    graph_path,
                    lib_path,
                    primitive_only,
                    S_bounds,
                    target_n,
                    maxNodes,
                    hc,
                    lc,
                    t,
                    output_path,
                ]
            )
            + "\n"
        )
