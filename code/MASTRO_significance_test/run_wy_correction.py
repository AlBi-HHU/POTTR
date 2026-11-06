import os
import argparse
import subprocess
import sys
import numpy as np
from tqdm import tqdm
from multiprocessing import Pool

# claude suggestion:
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))  # MASTRO_significance_test/
CODE_DIR = os.path.dirname(SCRIPT_DIR)  # Code/
POTTR_SCRIPT = os.path.join(CODE_DIR, 'run_POTTR.py')

QUANTILES = [0.01, 0.05, 0.1]


# claude suggestion:
def run_analysis(cmd):
    result = subprocess.run(cmd, cwd=CODE_DIR, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    return cmd, result.returncode, result.stderr


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", type=str, required=True, help="input file with graphs")
    parser.add_argument("-o", type=str, required=True, help="path to output")
    parser.add_argument("-m", type=int, required=True, help="number of resamples to create")
    parser.add_argument("-c", type=int, default=1, help="number of cores")
    parser.add_argument("-minp", type=str, required=True,
                        help="folder for the minimum pvalues; one subfolder k_<k> per k and a summary file")
    parser.add_argument("-par", type=int, default=1,
                        help="run the analysis in parallel: 0 = no, 1 = yes (def=1)")
    parser.add_argument("-kmin", type=int, default=2, help="smallest k to run (default 2)")
    parser.add_argument("-kmax", type=int, required=True, help="largest k to run")
    parser.add_argument("-pool", type=int, default=1, help="solution pool size as in the real POTTR run (default 1)")
    parser.add_argument("-rt", type=int, help="resolution threshold as in the real POTTR run (optional)")
    parser.add_argument("-rf", action='store_true', help="resolution frequency as in the real POTTR run (optional)")
    return parser


def create_resamples(args):
    # resamples do not depend on k and the seed is fixed, so create them only once
    if all(os.path.isdir(os.path.join(args.o, 'resample_' + str(i))) for i in range(args.m)):
        print('Resamples already exist, skip resampling')
        return
    cmd = [sys.executable, os.path.join(SCRIPT_DIR, 'create_random_permutations.py'),
           '-g', args.g, '-o', args.o, '-m', str(args.m), '-c', str(args.c)]
    print(' '.join(cmd))
    subprocess.run(cmd, check=True)


def main():
    args = get_parser().parse_args()

    # claude suggestion absolute paths
    # absolute paths: the POTTR runs are started from Code/, not from the current folder
    args.o = os.path.abspath(args.o)
    args.g = os.path.abspath(args.g)
    args.minp = os.path.abspath(args.minp)
    os.makedirs(args.o, exist_ok=True)

    create_resamples(args)

    # run_POTTR.py writes minp_values.txt into the given folder; one folder per k, cleared before the run
    minp_files = {}
    for k in range(args.kmin, args.kmax + 1):
        minp_dir = os.path.join(args.minp, 'k_' + str(k))
        os.makedirs(minp_dir, exist_ok=True)
        minp_files[k] = os.path.join(minp_dir, 'minp_values.txt')
        for f in [minp_files[k], minp_files[k] + '.lock']:
            if os.path.exists(f):
                os.remove(f)

    commands = []
    for k in range(args.kmin, args.kmax + 1):
        for i in range(args.m):
            path_resampled_db = os.path.join(args.o, 'resample_' + str(i))
            output_resample = os.path.join(args.o, 'pottr_output_k' + str(k), 'resample_' + str(i))
            command = [sys.executable, POTTR_SCRIPT,
                       '--dags', path_resampled_db,
                       '--output-path', output_resample,
                       '-k', str(k),
                       '-pool', str(args.pool),
                       '-c', '1',
                       '--minp', os.path.dirname(minp_files[k])]
            if args.rt is not None:
                command += ['-rt', str(args.rt)]
            if args.rf:
                command.append('-rf')
            commands.append(command)

    if args.par == 1:
        with Pool(args.c) as pool:
            results = list(tqdm(pool.imap_unordered(run_analysis, commands), total=len(commands)))
    else:
        results = [run_analysis(command) for command in tqdm(commands)]

    failures = [(command, code, err) for command, code, err in results if code != 0]
    n_failed = {k: 0 for k in range(args.kmin, args.kmax + 1)}
    for command, code, err in failures:
        k = int(command[command.index('-k') + 1])
        n_failed[k] += 1
        print('Failed:', ' '.join(command))
        print(err)

    # one line per k: number of min pvalues and the quantiles; significant if pval_POTTR < quantile
    summary_file = os.path.join(args.minp, 'minp_quantiles.csv')
    with open(summary_file, 'w') as f:
        f.write('k,n_minp,n_failed,' + ','.join('q_' + str(q) for q in QUANTILES) + '\n')
        for k in range(args.kmin, args.kmax + 1):
            minpvals = np.sort(np.loadtxt(minp_files[k], ndmin=1)) if os.path.exists(minp_files[k]) else np.array([])
            if len(minpvals) != args.m:
                print(f'Warning: k={k} has {len(minpvals)} min pvalues instead of {args.m}')
            quants = [str(minpvals[int(q * len(minpvals))]) if len(minpvals) > 0 else 'NA' for q in QUANTILES]
            f.write(f'{k},{len(minpvals)},{n_failed[k]},' + ','.join(quants) + '\n')
            print('k', k, 'quantiles', dict(zip(QUANTILES, quants)))

    print('Summary written to', summary_file)


if __name__ == "__main__":
    main()
