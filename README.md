OPTIMO
======

This repository contains supplementary software code for the book _Epicyclic Gearing: Optimization Techniques_ published by CRC Press.

## What is included? ##

- OPTIMO (python module) — Algorithms for Multi-Objective optimization.
- Problem sets for the case studies covered in the book.
- Sample problem definition demonstrating usage.

## Requirements ##

- Python 3.6 or newer
- Matplotlib (for plot generation)

## Problem sets ##

In order to run the included problem sets, execute `run.py` with the required options, or use `run.py --help` to see the complete list of options.

```sh
python run.py -a <algorithm> -p <problem>
```

Where,  
`<algorithm>` is one of `pso, sa, aco, nsga2`  
`<problem>` is one of `ch3, ch4, ch5, ch6a, ch6b, ch7a, ch7b`


## Usage ##

The following is the general structure of a problem definition based on the included `sample.py`.

Import the desired algorithm.

```py
from optimo.algorithm.pso import PSO
```

Provide the configuration information.

```py
config = {
  "num_vars": 2,
  "num_obj": 2,
  "num_constr": 2,
  "bounds_min": [0, 0],
  "bounds_max": [5, 3]
}
```

Provide the problem definition.

```py
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
```

Instantiate and run the solver, then save the results.

```py
alg = PSO(Problem, config)
alg.run()
alg.save()
```

## License ##

This software is released under the MIT License.
