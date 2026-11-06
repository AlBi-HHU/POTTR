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


def max_run(path):
    size = path.split('_')[-1]
    size = size.replace('/', '')
    if int(size) < 1024:
        return 5
    else:
        return 1


input_data = ['../../Data/simulated_dags_seed42/tracerx_128/',
              '../../Data/simulated_dags_seed42/tracerx_256/',
              '../../Data/simulated_dags_seed42/tracerx_512/',
              '../../Data/simulated_dags_seed42/tracerx_1024/',
              '../../Data/simulated_dags/tracerx_2048/'
]
output_path = '../../Output/benchmarking/new_benchmarking_runtime/'
os.makedirs(output_path, exist_ok=True)

result_file = output_path + 'average_runtimes_and_peak_mem.csv'
cpus = '128'

timelimit = 10 * 3600 # time limit of 10 hours in seconds

file_handle = open(result_file, 'w')
file_handle.write('Graph count,trajectory size,wall-clock time,user time,sys time,maximum resident set size\n')

for data in input_data:
    print(data)
    num_graphs = data.split('_')[-1].replace('/', '')
    cmd = ['/usr/bin/time', 'python', '../run_POTTR.py',
           '--output-path', output_path + num_graphs,
           '--dags', data,
           '-k', '2', '-c', cpus, '--parallelize']

    wall_times = []
    user_times = []
    sys_times = []
    mem_peaks = []
    timed_out = False

    runs = max_run(data)
    for i in range(runs):
        print('# Run', i, 'of', runs)
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
            print('Warning: Could not parse time output for run', i)

    if timed_out:
        print('Process {0}.{1} timed out after {2} seconds'.format(num_graphs, i, timelimit))
    else:
        wall_time =format(mean(wall_times), '.2f')
        user_time = format(mean(user_times), '.2f')
        sys_time = format(mean(sys_times), '.2f')
        mem_peak = mean(mem_peaks) / 1024  # convert KB to MB
        mem_peak = format(mem_peak, '.2f')
        file_handle.write(','.join([num_graphs, str(traj_size), wall_time, user_time, sys_time, mem_peak]) + '\n')
        file_handle.flush()

        if int(traj_size) == 0:
            break

file_handle.close()