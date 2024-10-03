Example of how to run continuity.py

```
>>> import RNA
>>> from continuity import ContinuousEvol
>>> from basic import Vienna
>>> c = ContinuousEvol(Vienna())
>>> c.evolution("((((....))))", RNA.random_string(12, "ACGU"), 10)
>>> c.population[0].history
[(0, 'GGGAAUGGUCCU', '((((....))))')]
```
