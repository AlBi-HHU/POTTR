import os
import re
import signal
import subprocess
from numpy import mean


def parse_wall_time(time_str):
    parts = time_str.split(':')
    if len(parts) == 2:    # m:ss.ss
        return float(parts[0]) * 60 + float(parts[1])
    elif len(parts) == 3:  # h:mm:ss
        return float(parts[0]) * 3600 + float(parts[1]) * 60 + float(parts[2])


input_data = ['../../Data/simulated_dags_seed42/tracerx_128/',
              '../../Data/simulated_dags_seed42/tracerx_256/'
]
output_path = '../../Output/benchmarks/new_benchmarking_k_runtime/'
os.makedirs(output_path, exist_ok=True)

result_file = output_path + 'k_average_runtimes_and_peak_mem_traj_sizes.csv'
cpus = '128'

timelimit = 10 * 3600 # time limit of 10 hours in seconds

file_handle = open(result_file, 'w')
file_handle.write('Graph count,k,trajectory size,wall-clock time,user time,sys time,maximum resident set size\n')

for data in input_data:
    print(data)
    num_graphs = data.split('_')[-1].replace('/', '')
    for k in range(2, int(num_graphs) + 1):
        print(k)
        cmd = ['/usr/bin/time', 'python', '../run_POTTR.py',
               '--output-path', output_path + num_graphs,
               '--dags', data,
               '-k', str(k), '-c', str(cpus), '--parallelize']

        wall_times = []
        user_times = []
        sys_times = []
        mem_peaks = []
        timed_out = False

        print('Start process')
        # start POTTR in a new process group
        process = subprocess.Popen(cmd,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE,
                                   start_new_session=True)

        try:
            stdout_, stderr_ = process.communicate(timeout=timelimit)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            process.communicate()
            timed_out = True
            break
        print('End process')

        # Decode the outputs
        stdout = stdout_.decode('utf-8')
        stderr = stderr_.decode('utf-8')

        traj_match = re.search(r'Size of a maximum trajectory for this instance is  (\d+)', stdout)
        wall_match = re.search(r'(\S+)elapsed', stderr)
        user_match = re.search(r'(\d+\.\d+)user', stderr)
        sys_match = re.search(r'(\d+\.\d+)system', stderr)
        mem_match = re.search(r'(\d+)maxresident', stderr)

        if wall_match and user_match and sys_match and mem_match:
            traj_size = traj_match.group(1)
            wall_times.append(parse_wall_time(wall_match.group(1)))
            user_times.append(float(user_match.group(1)))
            sys_times.append(float(sys_match.group(1)))
            mem_peaks.append(int(mem_match.group(1)))
        else:
            print(stderr)
            print('Warning: Could not parse time output for run', k)

        if timed_out:
            print('Process {0}.{1} timed out after {2} seconds'.format(num_graphs, k, timelimit))
        else:
            wall_time =format(mean(wall_times), '.2f')
            user_time = format(mean(user_times), '.2f')
            sys_time = format(mean(sys_times), '.2f')
            mem_peak = mean(mem_peaks) / 1024  # convert bytes to MB
            mem_peak = format(mem_peak, '.2f')
            file_handle.write(','.join([num_graphs, str(k), str(traj_size), wall_time, user_time, sys_time, mem_peak]) + '\n')
            file_handle.flush()

        if int(traj_size) == 0:
            break

file_handle.close()