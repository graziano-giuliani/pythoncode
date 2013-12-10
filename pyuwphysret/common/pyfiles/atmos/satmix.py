#!/usr/bin/env python
# satvap.py
import numpy as num
from satvap import satvap
from e2mr import e2mr

def satmix(p,T,Tconvert=None):
  """wsat = satmix(p,T,Tconvert)

compute saturation mixing ratio [g/kg] given reference pressure, 
p [mbar] and temperature, T [K].  If Tconvert input, the calculation uses 
the saturation vapor pressure over ice (opposed to over water) 
for temperatures less than Tconvert [K].

DCT, updated 3/5/00
  """
  # saturation pressure
  if (Tconvert is None):
    esat = satvap(T)
  else:
    esat = satvap(T,Tconvert)

  # saturation mixing ratio
  wsat = e2mr(p,esat)
  return wsat

if __name__ == '__main__':
  print(satmix.__doc__)
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
  t = t + 273.15
  w = satmix(p,t,253.15)
  print(w)
