#!/usr/bin/python

from optimo.algorithm.pso import PSO
from optimo.algorithm.sa import SA
from optimo.algorithm.aco import ACO
from optimo.algorithm.nsga2 import NSGA2
import argparse
import random

#############
# Arguments #
#############

problem_list = ['ch3', 'ch4', 'ch5', 'ch6a', 'ch6b', 'ch7a', 'ch7b']
algorithm_list = ['pso', 'sa', 'aco', 'nsga2']
format_list = ['png', 'pdf', 'svg']

parser = argparse.ArgumentParser()
parser.add_argument('-p', '--problem', choices=problem_list, type=str, required=True)
parser.add_argument('-a', '--algorithm', choices=algorithm_list, type=str, required=True)
parser.add_argument('-f', '--formats', nargs='+', choices=format_list, type=str, default=['png'])
parser.add_argument('-o', '--output', type=str, default='')
parser.add_argument('-d', '--debug', action='store_true')

args = parser.parse_args()

#######################
# Problem definition #
#######################

# Power Take-off Gearbox
if args.problem == 'ch3':
  from problemset.problem_ch3 import Problem, config

# Wind Turbine Generator Gearbox
elif args.problem == 'ch4':
  from problemset.problem_ch4 import Problem, config

# Helicopter Gearbox
elif args.problem == 'ch5':
  from problemset.problem_ch5 import Problem, config

# Electric Vehicle Gearbox - Speed 1
elif args.problem == 'ch6a':
  from problemset.problem_ch6a import Problem, config

# Electric Vehicle Gearbox - Speed 2
elif args.problem == 'ch6b':
  from problemset.problem_ch6b import Problem, config

# Wolfram Gearbox - Stage 1
elif args.problem == 'ch7a':
  from problemset.problem_ch7a import Problem, config

# Wolfram Gearbox - Stage 2
elif args.problem == 'ch7b':
  from problemset.problem_ch7b import Problem, config

#############
# Execution #
#############

random.seed(0.123)

if args.algorithm == 'pso':
  alg = PSO(Problem, config)
elif args.algorithm == 'sa':
  alg = SA(Problem, config)
elif args.algorithm == 'aco':
  alg = ACO(Problem, config)
elif args.algorithm == 'nsga2':
  alg = NSGA2(Problem, config)

if args.debug:
  alg.debug()

alg.run()

alg.set_path(args.output)
alg.set_formats(args.formats)
alg.save()

