#!/usr/bin/python

from optimo.algorithm import Algorithm
from optimo.utils import *
import random
import math
import copy

##############
# Parameters #
##############

class Particle:
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

class SA(Algorithm):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)

  def setup(self):
    self.population = [Particle() for i in range(self.num_pop)]
    self.archive = []

    self.step_size = 1.0
    self.temp = 100
    self.t = 0

    # Initialize population
    for p in self.population:
      self.init_particle(p)

    # Initialize archive
    self.update_archive()

  def init_particle(self, p):
    p.position = [random.uniform(self.bounds_min[i], self.bounds_max[i]) 
                  for i in range(self.num_vars)]

    self.evaluate(p.position, p.fitness, p.constraints)
    self.update_position(p)

    p.is_dominated = False
    p.crowding_dist = 1.0
    p.cof = 1.0

  def update_position(self, p):
    p.best_position = p.position[:]
    p.best_fitness = p.fitness[:]
    p.best_constraints = p.constraints[:]
    p.feasibility = Feasibility(p)

  def update_particle(self, p):
    for j in range(self.num_vars):
      r = random.uniform(-1.0, 1.0)
      p.position[j] = p.position[j] + r * self.step_size
      p.position[j] = Clamp(p.position[j], self.bounds_min[j], self.bounds_max[j])

    self.evaluate(p.position, p.fitness, p.constraints)

    if Feasibility(p) >= 0 and p.feasibility >= 0:
      if Dominates(p.fitness, p.best_fitness):
        self.update_position(p)
      else:
        t = self.temp * (1.0 - float(self.t) / self.num_iter)
        diff = [a - b for a,b in zip(p.fitness, p.best_fitness)]
        metropolis = math.exp(min(-sum(diff) / t, 0))
        if random.random() < metropolis:
          self.update_position(p)
    else:
      if Feasibility(p) > p.feasibility:
        self.update_position(p)

  def update_pop(self):
    for p in self.population:
      self.update_particle(p)

    self.t += 1

  def update_archive(self):
    self.archive = self.archive + copy.deepcopy(self.population)
    self.archive = FilterFeasible(self.archive)
    self.archive = FilterDistinct(self.archive)

    ComputeDominance(self.archive)
    self.archive = FilterNonDominated(self.archive)

    ComputeCrowdingDistance(self.archive)
    while len(self.archive) > self.num_archive:
      r = SelectMemberToDelete(self.archive)
      self.archive.pop(r)

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

