#!/usr/bin/python

from optimo.algorithm.pso import PSO

#######################
# Problem definition #
#######################

config = {
  "num_vars": 2,
  "num_obj": 2,
  "num_constr": 2,
  "bounds_min": [0, 0],
  "bounds_max": [5, 3]
}

# BNH (Binh and Korn)
def Problem(v, obj, con):
  x1, x2 = v[0], v[1]

  f1 = 4 * (x1 ** 2) + 4 * (x2 ** 2)
  f2 = ((x1 - 5) ** 2) + ((x2 - 5) ** 2)
  obj.append(f1)
  obj.append(f2)

  c1 = 25 - (((x1 - 5) ** 2) + (x2 ** 2))
  c2 = (((x1 - 8) ** 2) + ((x2 + 3) ** 2)) - 7.7
  con.append(c1)
  con.append(c2)

#############
# Execution #
#############

alg = PSO(Problem, config)
alg.run()
alg.save()

