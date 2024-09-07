#!/usr/bin/python

from optimo.algorithm import Algorithm
from optimo.utils import *
import random
import copy

##############
# Parameters #
##############

class Ant:
  path = []
  fitness = []
  constraints = []
  best_position = []
  best_fitness = []
  best_constraints = []
  feasibility = 0.0
  is_dominated = False
  crowding_dist = 1.0
  cof = 1.0

class Node:
  value = 0.0
  edges = []

class Edge:
  next_node = None
  pheromone = 1.0

##################
# Implementation #
##################

class ACO(Algorithm):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)

  def setup(self):
    self.population = [Ant() for i in range(self.num_pop)]
    self.archive = []

    self.divisions = 101
    self.rho = 0.1

    # Initialize nodes and edges
    self.init_nodes()
    self.init_edges()

    # Initialize ants
    for p in self.population:
      self.init_ant(p)

  def init_nodes(self):
    self.home = Node()
    self.home.value = 0.0
    self.home.edges = []

    self.nodes = []
    for i in range(self.num_vars):
      a, b = self.bounds_min[i], self.bounds_max[i]
      layer = []

      for j in range(self.divisions):
        n = Node()
        n.value = a + (b - a) * j / (self.divisions - 1)
        n.edges = []
        layer.append(n)

      self.nodes.append(layer)

  def init_edges(self):
    self.edges = []
    for j in range(self.divisions):
      e = Edge()
      e.next_node = self.nodes[0][j]
      e.pheromone = 1.0
      self.home.edges.append(e)
      self.edges.append(e)

    for i in range(self.num_vars-1):
      for n in self.nodes[i]:
        for j in range(self.divisions):
          e = Edge()
          e.next_node = self.nodes[i+1][j]
          e.pheromone = 1.0
          n.edges.append(e)
          self.edges.append(e)

  def init_ant(self, p):
    p.path = []
    p.fitness = []
    p.constraints = []
    p.best_position = []
    p.best_fitness = [1e8 for i in range(self.num_obj)]
    p.best_constraints = [-1e8 for i in range(self.num_constr)]
    p.feasibility = 0.0
    p.is_dominated = False
    p.crowding_dist = 1.0
    p.cof = 1.0

  def update_position(self, p):
    position = [e.next_node.value for e in p.path]
    p.best_fitness = p.fitness[:]
    p.best_constraints = p.constraints[:]
    p.best_position = position[:]
    p.feasibility = Feasibility(p)

  def update_ant(self, p):
    node = self.home
    p.path = []

    for j in range(self.num_vars):
      r = SelectEdge(node)
      edge = node.edges[r]
      node = edge.next_node
      p.path.append(edge)

    position = [e.next_node.value for e in p.path]
    self.evaluate(position, p.fitness, p.constraints)

    if Feasibility(p) >= 0 and p.feasibility >= 0:
      if Dominates(p.fitness, p.best_fitness):
        self.update_position(p)
      elif random.random() > 0.5:
        self.update_position(p)
    else:
      if Feasibility(p) > p.feasibility:
        self.update_position(p)

  def update_pheromone(self):
    for e in self.edges:
      e.pheromone = e.pheromone * (1.0 - self.rho)
      e.pheromone = max(e.pheromone, 1.0)

    for p in self.population:
      if Feasibility(p) >= 0:
        delta = 1.0 / sum(p.fitness)
        for e in p.path:
          e.pheromone += delta

    # Clear out path references here so deepcopy() takes a 
    # sane amount of time / space.
    for p in self.population:
      p.path = []

  def update_pop(self):
    for p in self.population:
      self.update_ant(p)

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
      self.update_pheromone()
      self.update_archive()

    self.finalize()

