#!/usr/bin/env python
import numpy
import supsmu

def SuperSmoother(x,y,w=None,iper=1,span=0.0,alpha=1.0):
  """smo = SuperSmoother(x,y,w,iper,span,alpha)

 Super-Smoother
 --------------

 Friedman J.H. (1984). A variable span smoother.
 Department of Statistics, Stanford University, Technical Report LCS5

 version 10/10/84.

 coded  and copyright (c) 1984 by:

   Jerome H. Friedman
   Department of Statistics
   and
   Stanford Linear Accelerator Center
   Stanford University
   all rights reserved.

Input:
  x = ordered abscissa values.
  y = corresponding ordinate (response) values.

Optional Input:
  w     : weight for each (x,y) observation.
  iper  : periodic variable flag.
           iper = 1 => x is ordered interval variable (default)
           iper = 2 => x is a periodic variable with values
                       in the range (0.0,1.0) and period 1.0
  span  : smoother span (fraction of observations in window)
           span=0.0 => automatic (variable) span selection (default)
  alpha : controles high frequency (small span) penality
          used with automatic span selection (bass tone control).
          alpha must be in the range (0.0:10.0)

Output:
  smo:  smoothed ordinate (response) values.

Note:
     for small samples (n < 40) or if there are substantial serial
     correlations between obserations close in x - value, then
     a prespecified fixed span smoother (span > 0) should be
     used. reasonable span values are 0.2 to 0.4.
  """
  iper = numpy.array(iper)
  span = numpy.array(span)
  alpha = numpy.array(alpha)
  smo = numpy.zeros(x.shape[0])
  if (w is None):
    w = numpy.ones(x.shape[0])
  supsmu.supsmu_mod(x,y,w,iper,span,alpha,smo)
  return smo

if (__name__ == '__main__'):
  import sys
  print(SuperSmoother.__doc__)
  x = numpy.linspace(-2,2,1001)
  y = numpy.sin(x)+numpy.random.randn(len(x))*0.1
  smo = SuperSmoother(x,y)

  if ( sys.version_info[0] == 2 ):
    from pylab import *
    plot(y)
    plot(smo)
    show()
  else:
    print(y)
    print(smo)
