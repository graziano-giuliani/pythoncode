#!/usr/bin/env python
# mr2dp.py
import numpy as num
from esice_goffgratch import esice_goffgratch
from eswat_goffgratch import eswat_goffgratch
from mr2e import mr2e

def mr2dp(p,t,mr,Tconvert=None):
  """dp = mr2dp(p,t,mr,Tconvert)

Compute dew point temperature given watervapor mass mixing ratio (g/kg).  
Uses Goff-Gratch formulation.

If input, the air temperature (t,K) is used to determine if
saturation over water (for t > Tconvert) or saturation over ice
(for t <= Tconvert) is used.

Two dew points are returned: dp1 is with RH defined as the ratio 
of water vapor partial pressure to saturation vapor pressure and
dp2 is with RH defined as the ratio of water vapor mixing ratio to 
saturation mixing ratio.

Notes: results not valid for dew points >= 370 K and  <= 160 K.

Reference: Goff-Gratch formulation from sixth revised 
      edition of Smithsonian Meteorology Tables.

DCT 3/6/01
  """
  # compute H2O partial pressure
  e = mr2e(p,mr)

  # interpolate (saturation pressure vs T) to desired pressure
  T = num.linspace(100.0,400.0,1501)
  esat_water = eswat_goffgratch(T)
  dp = num.interp(e,esat_water,T)

  # over ice for air temperature <= Tconvert
  if ( Tconvert is not None ):
    esat_ice = esice_goffgratch(T)
    ind = num.where(t <= Tconvert)
    dp[ind] = num.interp(e,esat_ice,T)[ind]
  return dp

if __name__ == '__main__':
  print(mr2dp.__doc__)
  p = num.array(
    ( 1012.0, 991.3, 969.1, 945.5, 920.4, 893.8, 865.7, 836.1, 805.1, 772.8,
       739.5, 705.2, 670.3, 635.0, 599.7, 564.5, 529.8, 495.7, 462.6, 430.7,
       400.0, 370.8, 343.0, 316.7, 292.0, 266.8, 247.2, 227.0, 208.2, 190.8,
       174.7, 159.9, 146.2, 133.6, 121.9, 111.3, 101.5,  92.6,  84.4,  76.9,
        70.0 ))
  t = num.array(
      ( 24.54, 23.16, 21.67, 20.23, 18.86, 17.49, 16.10, 14.69, 13.22, 11.52,
         9.53,  7.24,  4.80,  2.34,  0.04, -2.29, -4.84, -7.64,-10.66,-13.95,
       -17.54,-21.45,-25.58,-29.90,-34.33,-38.94,-43.78,-48.80,-53.94,-58.79,
       -63.27,-67.32,-70.74,-73.62,-75.74,-77.07,-77.43,-76.63,-75.06,-73.14,
       -71.43 ))
  t = t + 273.15
  r = num.array(
    ( 17.78, 16.92, 15.93, 14.87, 13.78, 12.70, 11.84, 10.96, 10.15,  9.31,
       8.46,  7.73,  7.05,  6.32,  5.62,  4.91,  4.10,  3.30,  2.67,  2.15,
       1.66,  1.26,  0.95,  0.68,  0.45,  0.28,  0.17,  0.10,  0.06,  0.04,
       0.02,  0.02,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,
       0.02 ))
  dp = mr2dp(p,t,r,253.15)
  print(dp)
