#!/usr/bin/env python
# tpwv.py
import numpy as num
from svpice import svpice
from svpwat import svpwat

def tpwv(p,r):
  """tpwvr = tpwv(p,r)
Calculates the total precipitable water vapor for a sounding

Inputs: (all vectors of same length)
    p = low-level layer average pressure(s) (mb)
    r = water vapor mixing ratio (kg/kg)

Outputs:
 tpwv = total precipitable water vapor (cm)

RLT, 010709
  """
  tpw = 0.0
  n = p.size
  pa = p*100.0    # convert from mb to Pa
  w = r*1000.0    # convert to g/kg
  for i in range(0,n-1):
    dp = pa[i] - pa[i+1]
    w_mean = 0.5 * (w[i] + w[i+1])
    tpw = tpw + w_mean * dp
  tpw = tpw/98.1 # convert from kg*m^-2 to cm
  return tpw

if __name__ == '__main__':
  print(tpwv.__doc__)
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
  r = r / 1000.0
  tpw = tpwv(p,r)
  print(tpw)
