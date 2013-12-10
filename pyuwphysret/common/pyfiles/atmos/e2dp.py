#!/usr/bin/env python
# dp2mr.py
import numpy as num
from eswat_goffgratch import eswat_goffgratch
from esice_goffgratch import esice_goffgratch

def e2dp(e,t,Tconvert=None):
  """dp = e2dp(e,t,Tconvert)

Compute dew point temperature given water vapor saturation 
pressure (e, mbar).  Uses Goff-Gratch formulation.

If input, the air temperature (t,K) is used to determine if
saturation over water (for t > Tconvert) or saturation over ice
(for t <= Tconvert) is used.

Notes: results not valid for dew points >= 370 K and  <= 160 K.

Reference: Goff-Gratch formulation from sixth revised 
      edition of Smithsonian Meteorology Tables.

DCT 3/6/01
  """
  T = num.linspace(100.0,400.0,1501)
  esat_water = eswat_goffgratch(T)
  dp = num.interp(e,esat_water,T)
  if ( Tconvert is not None ):
    esat_ice = esice_goffgratch(T)
    ind = num.where(t <= Tconvert)
    dp[ind] = num.interp(e,esat_ice,T)[ind]
  return dp

if __name__ == '__main__':
  from dp2e import dp2e
  print(e2dp.__doc__)
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
  dp = e2dp(e,t,253.15)
  print(dp)
