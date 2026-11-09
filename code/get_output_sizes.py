import os
import argparse
import networkx as nx
import csv
import time
from datetime import datetime


def process_file(file_path):
    graph_nodes = set()

    # Read the first line of the file
    with open(file_path, 'r') as f:
        line = f.readline().strip().split()

    # Process each part of the first line
    for e in line:
        e = e.replace('_', '')  # Handle special character replacement
        if '->-' in e:
            a, b = e.split('->-')
            graph_nodes.add(a)
            graph_nodes.add(b)
        elif '-?-' in e:
            a, b = e.split('-?-')
            graph_nodes.add(a)
            graph_nodes.add(b)
        elif '-/-' in e:
            a, b = e.split('-/-')
            graph_nodes.add(a)
            graph_nodes.add(b)

    graph_nodes.discard('0')

    # Return graph order excluding the root node
    return len(graph_nodes)


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-path', '-o', required=True, dest='path', type=str,
                        help='Path to store output files')
    parser.add_argument('--dags', '-d', required=True, dest='trees', type=str,
                        help='File or directory containing transitively closed DAGs')
    parser.add_argument('--max-k', '-k', required=True, dest='max_k', type=int)
    parser.add_argument('--pool-size', '-pool', dest='pool', type=int, default=1)
    parser.add_argument('--cores', '-c', dest='cores', type=int, default=128)
    parser.add_argument('--pottr-path', '-pottr', required=True, dest='pottr_path', type=str,
                        help='Path to directory containing POTTR code')

    return parser


def main():
    args = get_parser().parse_args()

    trees = args.trees
    out_dir = os.path.join(args.path, datetime.now().strftime("%Y-%m-%d_%H-%M-%S"))
    suffix = ''
    if args.pool:
        suffix = '_pool' + str(args.pool)
    cores = args.cores

    results = []
    k_counter = 2
    find_trajectory = True
    times = []
    while find_trajectory:
        out = os.path.join(out_dir, 'out_k' + str(k_counter) + suffix)
        os.makedirs(out, exist_ok=True)

        print(k_counter)
        cmd = ('python ' + args.pottr
               + '-o ' + out
               + ' --dags ' + trees
               + ' -k ' + str(k_counter)
               + ' -c ' + str(cores)
               + ' -parallel'
               + ' -pool ' + str(args.pool))
        start = time.time()
        os.system(cmd)
        end = time.time() - start
        times.append(('out_k' + str(k_counter), end))
        print(end)

        file_path = out + '/converted_graphs.txt'
        traj_size = process_file(file_path)
        results.append(('out_k' + str(k_counter), traj_size))

        if traj_size == 0 or k_counter > args.max_k:
            find_trajectory = False

        k_counter += 1

    # Sort results by folder name
    results.sort(key=lambda x: int(x[0].split('_k')[1]))

    # Write results to a CSV file
    output_file = os.path.join(out_dir, "graph_orders.csv")
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Folder", "Graph Order"])  # Header
        writer.writerows(results)

    output_times = os.path.join(out_dir, "times.csv")
    with open(output_times, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Folder", "Time [sec]"])  # Header
        writer.writerows(times)

    print(f"Results written to {output_file}")
    print(max(times))


if __name__ == "__main__":
    main()