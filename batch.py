#!/usr/bin/python

import os, sys

problem_list = ['ch3', 'ch4', 'ch5', 'ch6a', 'ch6b', 'ch7a', 'ch7b']
algorithm_list = ['pso', 'sa', 'aco', 'nsga2']
format_list = ['png', 'pdf', 'svg']

command = f'"{sys.executable}" run.py -f {' '.join(format_list)}'
os.makedirs('results', exist_ok=True)

for p in problem_list:
  for a in algorithm_list:
    path = os.path.join('results', p, a)
    os.makedirs(path, exist_ok=True)
    os.system(f'{command} -p {p} -a {a} -o "{path}" -d')
