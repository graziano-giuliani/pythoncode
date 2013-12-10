#!/usr/bin/env python
# rh2e.py
import numpy as num
from satvap import satvap
from satmix import satmix
from e2mr import e2mr

def rh2e(p,t,rh,Tconvert=None):
  """(e1,e2) = rh2e(p,t,rh,Tconvert)

determine H2O partial pressure (e,mbar) given
pressure (p,mbar), temperature (t,K), and
relative humidity (rh,%)

Two H2O partial pressures are returned: e1 is with RH defined as 
the ratio of water vapor partial pressure to saturation vapor
pressure and e2 is with RH defined as the ratio of water vapor 
mixing ratio to saturation mixing ratio.

if input, Tconvert is used as the temperature point to switch
from using saturation vapor pressure over water to over ice.

DCT 3/6/00
  """
  # ratio of water mass to dry air mass
  eps = 0.621970585

  # saturation pressure
  if (Tconvert is None):
    esat = satvap(t)   # Goff Gratch formulation, over water
    wsat = satmix(p,t) # Goff Gratch formulation, over water
  else:
    esat = satvap(t,Tconvert)   # Goff Gratch formulation, over water/ice
    wsat = satmix(p,t,Tconvert) # Goff Gratch formulation, over water/ice

  # H2O partial pressure w/ RH defined as ratio of pressures
  e1 = rh/100.0*esat
  # H2O partial pressure w/ RH defined as ratio of mixing ratios
  e2 = (rh/100.0*esat/(p-esat)*p)/(rh/100.0*esat/(p-esat)+1.0);
  return(e1,e2)

if __name__ == '__main__':
  from dp2e import dp2e
  from e2rh import e2rh
  print(rh2e.__doc__)
  t = num.array(
      ( 24.54, 23.16, 21.67, 20.23, 18.86, 17.49, 16.10, 14.69, 13.22, 11.52,
         9.53,  7.24,  4.80,  2.34,  0.04, -2.29, -4.84, -7.64,-10.66,-13.95,
       -17.54,-21.45,-25.58,-29.90,-34.33,-38.94,-43.78,-48.80,-53.94,-58.79,
       -63.27,-67.32,-70.74,-73.62,-75.74,-77.07,-77.43,-76.63,-75.06,-73.14,
       -71.43 ))
  t = t + 273.15
  td = num.array(
      ( 295.99569524, 294.88592297, 293.58149854, 292.11729779, 290.51490282,
        288.80633219, 287.25337561, 285.56579921, 283.86054795, 281.99074887,
        279.96863936, 278.00807838, 276.00353817, 273.74197577, 271.36371593,
        268.74827599, 265.5596088,  261.9472149,  258.46973102, 255.00425602,
        251.12242488, 247.15405877, 243.22262393, 238.86585074, 233.8823144,
        228.4539335,  223.20007008, 217.86176743, 212.95046128, 209.08799585,
        203.25047643, 202.6535621,  197.18886555, 196.61856765, 196.0340168,
        195.44634221, 194.83729251, 194.21361376, 193.57543455, 192.93607596,
        196.90293301))
  p = num.array(
      ( 1012.0, 991.3, 969.1, 945.5, 920.4, 893.8, 865.7, 836.1, 805.1, 772.8,
         739.5, 705.2, 670.3, 635.0, 599.7, 564.5, 529.8, 495.7, 462.6, 430.7,
         400.0, 370.8, 343.0, 316.7, 292.0, 266.8, 247.2, 227.0, 208.2, 190.8,
         174.7, 159.9, 146.2, 133.6, 121.9, 111.3, 101.5,  92.6,  84.4,  76.9,
          70.0 ))
  e = dp2e(t,td,253.15)
  (rh1,rh2) = e2rh(p,t,e,253.15)
  (e1,e2) = rh2e(p,t,rh1,253.15)
  print(e1)
  (e1,e2) = rh2e(p,t,rh2,253.15)
  print(e2)
