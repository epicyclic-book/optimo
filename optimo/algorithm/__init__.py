#!/usr/bin/python

from optimo.utils import *
from optimo.report import *
import os
import math
import random
import copy

##################
# Algorithm Base #
##################

class Algorithm:
  
  def __init__(self, problem, config):
    super().__init__()

    self.problem = problem

    self.num_vars = config.get('num_vars', 0)
    self.num_obj = config.get('num_obj', 0)
    self.num_constr = config.get('num_constr', 0)
    self.bounds_min = config.get('bounds_min', [])
    self.bounds_max = config.get('bounds_max', [])
    self.ref_data = config.get('ref_data', [])

    self.num_iter = config.get('num_iter', 100)
    self.num_pop = config.get('num_pop', 200)
    self.num_archive = config.get('num_archive', 200)

    self.obj_names = config.get('obj_names', [])
    self.constr_names = config.get('constr_names', [])
    self.var_names = config.get('var_names', [])
    self.var_types = config.get('var_types', [])
    self.plot_data = config.get('plot_data', [])

    self.prepare_config()

    self.formats = ['png']
    self.path = ''

    self.population = []
    self.archive = []

  def prepare_config(self):
    if not self.obj_names:
      self.obj_names = [f'f{i+1}' for i in range(self.num_obj)]
    if not self.constr_names:
      self.constr_names = [f'c{i+1}' for i in range(self.num_constr)]
    if not self.var_names:
      self.var_names = [f'x{i+1}' for i in range(self.num_vars)]
    if not self.var_types:
      self.var_types = ['float' for i in range(self.num_vars)]
    if not self.plot_data:
      self.plot_data = [self.obj_names[:3]]
    if not self.ref_data:
      self.ref_data = self.bounds_min[:]

  def set_path(self, path):
    if os.path.isdir(path):
      self.path = path

  def set_formats(self, formats):
    format_list = {'png', 'pdf', 'svg'}
    self.formats = set(formats) & format_list

  def evaluate(self, v, obj, con):
    obj.clear()
    con.clear()

    try:
      self.problem(v, obj, con)
    except (ValueError, ArithmeticError):
      obj.clear()
      con.clear()

      obj.extend([ 1e8 for i in range(self.num_obj)])
      con.extend([-1e8 for i in range(self.num_constr)])

  def debug(self):
    v = self.ref_data
    obj, con = [], []
    self.evaluate(v, obj, con)

    print('Reference data:')
    print('Variables:', v)
    print('Objectives:', obj)
    print('Constraints:', con)

    assert len(self.bounds_min) == self.num_vars
    assert len(self.bounds_max) == self.num_vars
    assert len(self.obj_names) == self.num_obj
    assert len(self.constr_names) == self.num_constr
    assert len(self.var_names) == self.num_vars
    assert len(self.var_names) == len(self.var_types)
    assert all([obj in self.obj_names for p in self.plot_data for obj in p])
    assert all([2 <= len(p) <= 3 for p in self.plot_data])

    assert len(self.var_names) == len(v)
    assert len(self.obj_names) == len(obj)
    assert len(self.constr_names) == len(con)

  def run(self):
    pass

  def save(self):
    self.report()
    self.plot()

  def report(self):
    pop = self.archive
    header = self.obj_names + self.constr_names + self.var_names
    output = os.path.join(self.path, 'results')
    Report(pop, header, self.var_types, output)

  def plot(self):
    pop = self.archive
    data = copy.deepcopy(pop)

    for n, header in enumerate(self.plot_data):
      indices = [self.obj_names.index(p) for p in header]

      for d, p in zip(data, pop):
        d.best_fitness = [p.best_fitness[i] for i in indices]

      ComputeDominance(data)
      dataset = (
        (FilterDominated(data), 'blue', 'o'),
        (FilterNonDominated(data), 'red', '^')
      )

      output = os.path.join(self.path, f'plot{n+1}')
      Plot(dataset, header, output, self.formats)

#####################
# Genetic Algorithm #
#####################

class GeneticAlgorithm(Algorithm):

  def __init__(self, *args, **kwargs):
    super().__init__(*args, **kwargs)
    self.num_pop -= self.num_pop % 4
    self.num_archive = self.num_pop

  # Tournament selection
  def selection(self):
    self.population = []

    pop = self.archive
    indices = [i for i in range(self.num_pop)]

    for i in range(2):
      random.shuffle(indices)
      for i in range(0, self.num_pop, 4):
        p1 = self.tournament(pop[indices[i  ]], pop[indices[i+1]])
        p2 = self.tournament(pop[indices[i+2]], pop[indices[i+3]])
        q1, q2 = self.crossover(p1, p2)

        self.population.append(q1)
        self.population.append(q2)

  def tournament(self, p, q):
    if p.feasibility >= 0 and q.feasibility >= 0:
      if Dominates(p.best_fitness, q.best_fitness):
        return p
      if Dominates(q.best_fitness, p.best_fitness):
        return q
    else:
      if p.feasibility > q.feasibility:
        return p
      if q.feasibility > p.feasibility:
        return q

    if p.crowding_dist > q.crowding_dist:
      return p
    if q.crowding_dist > p.crowding_dist:
      return q

    if random.random() > 0.5:
      return p
    else:
      return q

  # Simulated binary crossover (K. Deb et al. 2007)
  # https://doi.org/10.1145/1276958.1277190
  def crossover(self, p1, p2, prob=0.8, eta=10):
    q1 = copy.deepcopy(p1)
    q2 = copy.deepcopy(p2)

    if random.random() > prob:
      return q1, q2

    for i in range(self.num_vars):
      if random.random() > 0.5:
        continue

      y1 = p1.position[i]
      y2 = p2.position[i]
      if abs(y1 - y2) < 1e-8:
        continue

      c = [y1, y2]
      y1, y2 = sorted(c)
      yl, yu = self.bounds_min[i], self.bounds_max[i]

      rand = random.random()
      exp = 1.0 / (eta + 1.0)
      delta = y2 - y1

      # Child 1
      beta = 1.0 + 2.0 * (y1 - yl) / delta
      alpha = 2.0 - math.pow(beta, -(eta + 1.0))

      alpha *= rand
      if alpha <= 1.0:
        betaq = math.pow(alpha, exp)
      else:
        betaq = math.pow(1.0 / (2.0 - alpha), exp);

      c[0] = 0.5 * ((y1 + y2) - betaq * delta)

      # Child 2
      beta = 1.0 + 2.0 * (yu - y2) / delta
      alpha = 2.0 - math.pow(beta, -(eta + 1.0))

      alpha *= rand
      if alpha <= 1.0:
        betaq = math.pow(alpha, exp)
      else:
        betaq = math.pow(1.0 / (2.0 - alpha), exp);

      c[1] = 0.5 * ((y1 + y2) + betaq * delta)
      
      random.shuffle(c)
      q1.position[i] = Clamp(c[0], yl, yu)
      q2.position[i] = Clamp(c[1], yl, yu)

    return q1, q2

  # Polynomial mutation (K. Deb et al. 2007)
  # https://doi.org/10.1145/1276958.1277190
  def mutation(self, p, prob=0.5, eta=20):
    for i in range(self.num_vars):
      if random.random() > prob:
        continue

      y = p.position[i]
      yl, yu = self.bounds_min[i], self.bounds_max[i]
      if abs(yu - yl) < 1e-8:
        continue

      delta1 = (y - yl) / (yu - yl)
      delta2 = (yu - y) / (yu - yl)
      rand = random.random()
      exp = 1.0 / (eta + 1.0)

      if rand <= 0.5:
        x = math.pow(1.0 - delta1, eta + 1.0)
        x = 2.0 * rand + (1.0 - 2.0 * rand) * x
        deltaq = math.pow(x, exp) - 1.0
      else:
        rand = 1.0 - rand
        x = math.pow(1.0 - delta2, eta + 1.0)
        x = 2.0 * rand + (1.0 - 2.0 * rand) * x
        deltaq = 1.0 - math.pow(x, exp)

      c = y + deltaq * (yu - yl)
      p.position[i] = Clamp(c, yl, yu)

