#!/usr/bin/env python
# parcel_thetae.py
import numpy as num
from wvmr import wvmr
from temp_at_lcl import temp_at_lcl

def parcel_thetae(p,t,td,p0):
  """thetaE_parcel = parcel_thetae(p,t,td,p0)
Calculates the parcel equivalent potential temperature

Inputs: (all vectors of same length)
    p = low-level layer average pressure(s) (mb)
    t = dry bulb temperature(s) (K)
   td = dewpoint temperature(s) (K)
   p0 = surface pressure(s) (mb)

Outputs:
 thetae_parcel = parcel equivalent potential temperature (K)

RLT, 010709
  """
  theta = t*pow(p0/p,0.28541)
  thetaE_parcel = theta*num.exp(2481.9e-3*wvmr(p,t,td)/temp_at_lcl(t,td))
  return thetaE_parcel

if __name__ == '__main__':
  from tdew import tdew
  print(parcel_thetae.__doc__)
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
  p0 = 1013.0
  thetae = parcel_thetae(p,t,td,p0)
  print(thetae)
