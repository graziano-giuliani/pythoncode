#!/usr/bin/env python
# Lu.py
from scipy.linalg import lu

def Lu(A):
  """(L,U,P) = Lu(A)
  Compute pivoted LU decompostion of a matrix.

  RETURNS:
             L,U,P

      where A = PLU
  """
  (P,L,U) = lu(A)
  return ( L,U,P )

if ( __name__ == '__main__'):
  import numpy
  print(Lu.__doc__)
  A = numpy.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
  (L,U,P) = Lu(A)
  print(L)
  print(U)
  print(P)
  print(sum(sum(numpy.dot(P,A)-numpy.dot(L,U))))
