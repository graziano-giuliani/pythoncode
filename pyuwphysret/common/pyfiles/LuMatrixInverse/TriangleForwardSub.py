#!/usr/bin/env python
# TriangleForwardSub.py
from scipy.linalg import solve

def TriangleForwardSub(L,b):
  """C = TriangleForwardSub(L,b)
  Solve linear system  LC = b
  """
  C = solve(L,b)
  return C

if ( __name__ == '__main__'):
  import numpy
  print(TriangleForwardSub.__doc__)
  b = numpy.array([1.0,2.0,3.0])
  L = numpy.array(([1.0,0.0,0.0],[2.0,1.0,0.0],[3.0,4.0,1.0]))
  C = TriangleForwardSub(L,b)
  print(C)
