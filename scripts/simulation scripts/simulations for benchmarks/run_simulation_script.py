import os

start = 128
max = 3000
num = start

while num <= max:
    call = ('python simulated_trees_real_distribution.py '
            '-i ../../Data/data_tracerx/all_trees_hugo+type_drivers_only '
            '-o ../../Data/simulated_dags_seed42/tracerx '
            '-n ' + str(num) + ' '
            '-t 5')
    print(call)
    os.system(call)
    num *= 2