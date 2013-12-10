#!/usr/bin/env python
# svpice.py
import numpy as num

def svpice(t):
  """e = svpice(t)
Calculates the water vapor mixing ratio

Inputs: (all vectors of same length)
    t = dry bulb temperature(s) (K)

Outputs:
    e = saturation vapor pressure with respect to a plane surface of ice (mb)

RLT, 010709
  """
  A0=0.7859063157e0
  A1=0.357924232e-1
  A2=-0.1292820828e-3
  A3=0.5937519208e-6
  A4=0.4482949133e-9
  A5=0.2176664827e-10
  T = t - 273.16
  e = pow(10.0,A0+T*(A1 + T*(A2 + T*(A3 + T*(A4 + T*A5)))))
  return e

if __name__ == '__main__':
  print(svpice.__doc__)
  t = num.array(
      ( 24.54, 23.16, 21.67, 20.23, 18.86, 17.49, 16.10, 14.69, 13.22, 11.52,
         9.53,  7.24,  4.80,  2.34,  0.04, -2.29, -4.84, -7.64,-10.66,-13.95,
       -17.54,-21.45,-25.58,-29.90,-34.33,-38.94,-43.78,-48.80,-53.94,-58.79,
       -63.27,-67.32,-70.74,-73.62,-75.74,-77.07,-77.43,-76.63,-75.06,-73.14,
       -71.43 ))
  t = t + 273.15
  e = svpice(t)
  print(t)
  print(e)
