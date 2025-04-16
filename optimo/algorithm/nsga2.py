#!/usr/bin/python

from optimo.algorithm import GeneticAlgorithm
from optimo.utils import *
import random
import copy

##############
# Parameters #
##############

class Individual:
  position = []
  fitness = []
  constraints = []
  best_position = []
  best_fitness = []
  best_constraints = []
  feasibility = 0.0
  is_dominated = False
  crowding_dist = 1.0
  cof = 1.0

##################
# Implementation #
##################

class NSGA2(GeneticAlgorithm):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)

  def setup(self):
    self.population = [Individual() for i in range(self.num_pop)]
    self.archive = []

    # Initialize population
    for p in self.population:
      self.init_individual(p)

    # Initialize archive
    self.update_archive()

  def init_individual(self, p):
    p.position = [random.uniform(self.bounds_min[i], self.bounds_max[i]) 
                  for i in range(self.num_vars)]

    self.update_individual(p)

    p.is_dominated = False
    p.crowding_dist = 1.0
    p.cof = 1.0

  def update_position(self, p):
    p.best_position = p.position[:]
    p.best_fitness = p.fitness[:]
    p.best_constraints = p.constraints[:]
    p.feasibility = Feasibility(p)
    
  def update_individual(self, p):
    self.evaluate(p.position, p.fitness, p.constraints)
    self.update_position(p)

  def update_pop(self):
    self.selection()

    for p in self.population:
      self.mutation(p)
      self.update_individual(p)

  def update_archive(self):
    self.archive = self.archive + copy.deepcopy(self.population)
    self.archive = self.fill_archive(self.archive)

  # Naive extraction of non-dominated fronts
  def get_fronts(self, pop):
    while len(pop) > 0:
      ComputeDominance(pop)
      front = FilterNonDominated(pop)

      ComputeCrowdingDistance(front)
      pop = FilterDominated(pop)
      yield front

  # Fast non-dominated sorting (K. Deb et al. 2002)
  # https://doi.org/10.1109/4235.996017
  def get_fronts_fnds(self, pop):
    size = len(pop)
    dominance_set = [[] for i in range(size)]
    dominance_count = [0 for i in range(size)]

    def store(i, j):
      dominance_set[i].append(j)
      dominance_count[j] += 1

    ComputeDominance(pop, store)

    current_front = []
    for i in range(size):
      if dominance_count[i] == 0:
        current_front.append(i)

    while len(current_front) > 0:
      front = [pop[i] for i in current_front]
      ComputeCrowdingDistance(front)
      yield front

      next_front = []
      for i in current_front:
        for j in dominance_set[i]:
          dominance_count[j] -= 1
          if dominance_count[j] == 0:
            next_front.append(j)

      current_front = next_front

  def fill_archive(self, pop):
    archive = []
    for front in self.get_fronts_fnds(pop):
      if len(archive) + len(front) <= self.num_pop:
        archive.extend(front)
      else:
        front.sort(reverse=True, key=lambda p: p.crowding_dist)
        size = self.num_pop - len(archive)
        archive.extend(front[:size])
        break

    return archive

  def finalize(self):
    ComputeCOF(self.archive)
    self.archive.sort(key=lambda p: p.cof)

  def run(self):
    self.setup()

    for i in range(self.num_iter):
      print(f'Iteration: {i+1}')
      self.update_pop()
      self.update_archive()

    self.finalize()

