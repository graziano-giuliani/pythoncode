#!/usr/bin/env python
# MatrixInverse.py
from scipy.linalg import inv

def MatrixInverse(A):
  """Ainv = MatrixInverse(A)
  Compute the inverse of a matrix.

  RETURNS:
             Ainv

      where Ainv*A = I
  """
  Ainv = inv(A)
  return Ainv

if ( __name__ == '__main__'):
  import numpy
  print(MatrixInverse.__doc__)
  A = numpy.array([[1.0, 2.0], [3.0, 4.0]])
  Ainv = MatrixInverse(A)
  print(Ainv)
  print(numpy.dot(A,Ainv))
