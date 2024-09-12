#!/usr/bin/python

from matplotlib import pyplot
import math
import csv

##############
# Parameters #
##############

font = {
  'family': 'serif',
  'color': 'darkred',
  'weight': 'bold',
  'size': 12
}

######################
# Plotting functions #
######################

def Create2DPlot(header):
  fig = pyplot.figure()
  ax = fig.add_subplot()
  ax.set_xlabel(header[0], fontdict=font)
  ax.set_ylabel(header[1], fontdict=font)
  return fig, ax

def Create3DPlot(header):
  fig = pyplot.figure()
  ax = fig.add_subplot(projection='3d')
  ax.set_xlabel(header[0], fontdict=font)
  ax.set_ylabel(header[1], fontdict=font)
  ax.set_zlabel(header[2], fontdict=font)
  return fig, ax

def Plot(dataset, header, output, formats=['png']):
  n = len(header)
  if not 2 <= n <= 3:
    return

  if n == 2:
    fig, ax = Create2DPlot(header)
  else:
    fig, ax = Create3DPlot(header)

  dataset = [s for s in dataset if s[0]]
  for data, color, marker in dataset:
    points = zip(*(map(abs, p.best_fitness[:n]) for p in data))
    ax.scatter(*points, c=color, marker=marker)

  for fmt in formats:
    fig.savefig(f'{output}.{fmt}', dpi=300)

####################
# Report solutions #
####################

def Report(pop, header, types, output):
  with open(f'{output}.csv', 'w', newline='') as f:
    fieldnames = tuple(header) + ('Feasibility', 'Crowding Distance', 'COF')
    writer = csv.DictWriter(f, fieldnames=fieldnames, quoting=csv.QUOTE_ALL)
    writer.writeheader()

    for p in pop:
      process = lambda v, t: math.ceil(v) if t == 'int' else v
      position = list(map(process, p.best_position, types))

      fitness = list(map(abs, p.best_fitness))
      feasibility = int(p.feasibility >= 0)

      values = fitness + p.best_constraints + position
      values += [feasibility, p.crowding_dist, p.cof]

      row = dict(zip(fieldnames, values))
      writer.writerow(row)
