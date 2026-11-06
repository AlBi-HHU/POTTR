import os
import subprocess
import sys

TOY_DIR = os.path.dirname(os.path.abspath(__file__))      # Code/toy_example
CODE_DIR = os.path.dirname(TOY_DIR)                       # Code

x = '4'
input_dir = os.path.join(TOY_DIR, 'toy_graphs_resolved/' + x + '_resolved_edges/')
output_path = os.path.join(TOY_DIR, 'toy_graphs_resolved/' + x + '_resolved_edges/output/')
run_pottr = os.path.join(CODE_DIR, 'run_POTTR.py')

files = [f for f in os.listdir(input_dir) if f.endswith('.txt')]
sorted_files = sorted(files, key=lambda x: int(x.split('_')[0]))

for file in sorted_files:
  x = file.split('_')[0]
  cmd = [sys.executable, run_pottr,
         '--output-path', os.path.join(output_path, x),
         '--dags', os.path.join(input_dir, file),
         '-k', '2', '-c', '8', '--parallelize', '-v']
  print(' '.join(cmd))
  subprocess.run(cmd, check=True, cwd=CODE_DIR)