#!/usr/bin/env python
# TriangleBackwardSub.py
from scipy.linalg import solve

def TriangleBackwardSub(U,b):
  """C = TriangleBackwardSub(U,b)
  Solve linear system  UC = b
  """
  C = solve(U,b)
  return C

if ( __name__ == '__main__'):
  import numpy
  print(TriangleBackwardSub.__doc__)
  b = numpy.array([1.0,2.0,3.0])
  U = numpy.array(([2.0,2.0,1.0],[0.0,1.0,4.0],[0.0,0.0,3.0]))
  C = TriangleBackwardSub(U,b)
  print(C)
