#!/usr/bin/env python
# tdew.py
import numpy as num
from rh import rh

def tdew(p,t,r):
  """td = tdew(p,t,r)
Calculates dewpoint temperature of a sounding using the Arden Buck equation

Inputs: (arrays)
    p  = atmospheric pressure (mb) at N levels
    t  = dry bulb temperatures (K) at N levels
    r  = water vapor mixing ratio (kg H2O per kg dry air) at N levels

Outputs:
    td = Dewpoint temperature (K)
"""
  t0 = 273.15     # C-to-K or K-to-C linear conversion (deg)
  a = 6.112      # mb
  b = 18.678
  c = 257.14     # Celsius
  d = 234.5      # Celsius

  # water vapor pressure ev and saturation vapor pressure es
  tc = t - t0           # temperature in celsius
  rhum = rh(p,t,r)
  gamma_m = num.log(rhum*num.exp((b-tc/d)*(tc/(c+tc))))
  td = (c*gamma_m)/(b-gamma_m)
  return td+t0

if __name__ == '__main__':
  print(tdew.__doc__)
  t = num.array(
      ( 24.54, 23.16, 21.67, 20.23, 18.86, 17.49, 16.10, 14.69, 13.22, 11.52,
         9.53,  7.24,  4.80,  2.34,  0.04, -2.29, -4.84, -7.64,-10.66,-13.95,
       -17.54,-21.45,-25.58,-29.90,-34.33,-38.94,-43.78,-48.80,-53.94,-58.79,
       -63.27,-67.32,-70.74,-73.62,-75.74,-77.07,-77.43,-76.63,-75.06,-73.14,
       -71.43 ))
  p = num.array(
      ( 1012.0, 991.3, 969.1, 945.5, 920.4, 893.8, 865.7, 836.1, 805.1, 772.8,
         739.5, 705.2, 670.3, 635.0, 599.7, 564.5, 529.8, 495.7, 462.6, 430.7,
         400.0, 370.8, 343.0, 316.7, 292.0, 266.8, 247.2, 227.0, 208.2, 190.8,
         174.7, 159.9, 146.2, 133.6, 121.9, 111.3, 101.5,  92.6,  84.4,  76.9,
          70.0 ))
  r = num.array(
      ( 17.78, 16.92, 15.93, 14.87, 13.78, 12.70, 11.84, 10.96, 10.15,  9.31,
         8.46,  7.73,  7.05,  6.32,  5.62,  4.91,  4.10,  3.30,  2.67,  2.15,
         1.66,  1.26,  0.95,  0.68,  0.45,  0.28,  0.17,  0.10,  0.06,  0.04,
         0.02,  0.02,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,
         0.02 ))
  t = t + 273.15
  r = r / 1000.0
  td = tdew(p,t,r)
  print(td)
