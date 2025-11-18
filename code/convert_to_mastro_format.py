import argparse
import pandas as pd


def get_parser():
    """Get parser object for combinatorial vaccine design."""
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', '-i', required=True, dest='input', type=str,
                        help='Input counts and edges to convert to mastro format')
    parser.add_argument('--output_dir', '-o', required=True, dest='dir', type=str,
                        help='Output directory to store max trajectories and support')
    return parser


def convert(input_file, dir):
    df = pd.read_csv(input_file)
    output = dir + '/converted_graphs.txt'
    with open(output, 'w') as file:
        for row in df[['Support', 'Edges', 'Supporting Graphs']].iterrows():
            support = row[1]['Support']
            edges = row[1]['Edges']
            supporting_graphs = row[1]['Supporting Graphs']

            file.write(edges + ' (' + str(support) + ')\n')
            file.write(supporting_graphs.replace('L', '') + '\n')
    return output


if __name__ == '__main__':
    args = get_parser().parse_args()
    convert(args.input, dir)
