#!/usr/bin/python

import random
import itertools

####################
# Helper Functions #
####################

def Clamp(value, lower, upper):
  return max(lower, min(value, upper))

def RouletteWheelSelection(weights):
  total = sum(weights)
  weights = [w / total for w in weights]
  weights = itertools.accumulate(weights)

  r = random.random()
  less = [int(w < r) for w in weights]
  return sum(less)

def SelectGlobalOptimum(pop):
  weights = [p.crowding_dist for p in pop]
  return RouletteWheelSelection(weights)

def SelectMemberToDelete(pop):
  weights = [1.0 / p.crowding_dist for p in pop]
  return RouletteWheelSelection(weights)

def SelectEdge(node):
  weights = [e.pheromone for e in node.edges]
  return RouletteWheelSelection(weights)

def Dominates(p, q):
  less_eq = [i <= j for i,j in zip(p,q)]
  less = [i < j for i,j in zip(p,q)]
  return all(less_eq) and any(less)

def Feasibility(p):
  return sum([c for c in p.constraints if c < 0])

def ComputeDominance(pop, store = None):
  for p in pop:
    p.is_dominated = False

  size = len(pop)
  for i in range(size):
    for j in range(i+1, size):
      p, q = pop[i], pop[j]
      if p.feasibility >= 0 and q.feasibility >= 0:
        if Dominates(p.best_fitness, q.best_fitness):
          q.is_dominated = True
          store and store(i, j)
        if Dominates(q.best_fitness, p.best_fitness):
          p.is_dominated = True
          store and store(j, i)
      else:
        if p.feasibility > q.feasibility:
          q.is_dominated = True
          store and store(i, j)
        if q.feasibility > p.feasibility:
          p.is_dominated = True
          store and store(j, i)

def ComputeCrowdingDistance(pop):
  for p in pop:
    p.crowding_dist = 1.0

  size = len(pop)
  if size <= 2:
    return

  nobj = len(pop[0].best_fitness)
  if nobj < 1:
    return

  for p in pop:
    p.crowding_dist = 0.0

  for i in range(nobj):
    pop.sort(key=lambda p: p.best_fitness[i])
    delta = pop[-1].best_fitness[i] - pop[0].best_fitness[i]
    if abs(delta) < 1e-8:
      continue

    for j in range(size):
      a = Clamp(j-1, 0, size-1)
      b = Clamp(j+1, 0, size-1)
      p, a, b = pop[j], pop[a], pop[b]
      p.crowding_dist += (b.best_fitness[i] - a.best_fitness[i]) / delta

  for p in pop:
    p.crowding_dist /= nobj
    p.crowding_dist = Clamp(p.crowding_dist, 1e-8, 1.0)

def ComputeCOF(pop):
  for p in pop:
    p.cof = 0.0

  size = len(pop)
  if size < 2:
    return

  nobj = len(pop[0].best_fitness)
  if nobj < 1:
    return

  points = list(zip(*(p.best_fitness for p in pop)))
  cof_min = [min(p) for p in points]
  cof_max = [max(p) for p in points]

  for p in pop:
    for f, a, b in zip(p.best_fitness, cof_min, cof_max):
      if abs(b - a) >= 1e-8:
        p.cof += (f - a) / (b - a)

    p.cof /= nobj

def FilterFeasible(pop):
  return [p for p in pop if p.feasibility >= 0]

def FilterDominated(pop):
  return [p for p in pop if p.is_dominated]

def FilterNonDominated(pop):
  return [p for p in pop if not p.is_dominated]

def FilterDistinct(pop):
  out_pop = []
  added = set()

  for p in pop:
    key = hash(tuple(p.best_fitness))
    if key not in added:
      out_pop.append(p)
      added.add(key)

  return out_pop
