#!/usr/bin/env python
# svpwat.py
import numpy as num

def svpwat(t):
  """e = svpwat(t)
Calculates the water vapor mixing ratio

Inputs: (all vectors of same length)
    t = dry bulb temperature(s) (K)

Outputs:
    e = saturation vapor pressure with respect to a plane surface of ice (mb)

RLT, 010710
  """
  A0 = 0.999996876e0
  A1 = -0.9082695004e-2
  A2 = 0.7873616869e-4
  A3 = -0.6111795727e-6
  A4 = 0.4388418740e-8
  A5 = -0.2988388486e-10
  A6 = 0.2187442495e-12
  A7 = -0.1789232111e-14
  A8 = 0.1111201803e-16
  A9 = -0.3099457145e-19
  B = 0.61078e+1
  T = t - 273.16
  E = A0 + T*(A1 + T*(A2 + T*(A3 + T*(A4 + T*(A5 + T*(A6 + T*(A7 +
           T*(A8 + T*A9))))))))
  E = B/pow(E,8)
  return E

if __name__ == '__main__':
  print(svpwat.__doc__)
  t = num.array(
      ( 24.54, 23.16, 21.67, 20.23, 18.86, 17.49, 16.10, 14.69, 13.22, 11.52,
         9.53,  7.24,  4.80,  2.34,  0.04, -2.29, -4.84, -7.64,-10.66,-13.95,
       -17.54,-21.45,-25.58,-29.90,-34.33,-38.94,-43.78,-48.80,-53.94,-58.79,
       -63.27,-67.32,-70.74,-73.62,-75.74,-77.07,-77.43,-76.63,-75.06,-73.14,
       -71.43 ))
  t = t + 273.15
  e = svpwat(t)
  print(t)
  print(e)
