#!/usr/bin/env python
# dp2e.py
import numpy as num
from esice_goffgratch import esice_goffgratch
from eswat_goffgratch import eswat_goffgratch

def dp2e(t,dp,Tconvert=None):
  """e = dp2e(t,dp,Tconvert)

function e = dp2e(t,dp,Tconvert)

compute H2O partial pressure (e,mabr) given total pressure p (mb), 
air temperature t (K), and dew point temperature (K).

if input, Tconvert is used as the AIR temperature to switch
from using saturation vapor pressure over water to over ice.

dct 3/6/2000
  """
  # vapor pressures computed at dew point temperature over water and ice
  # initialize water vapor partial pressure as over water for all
  # air temperatures
  e = eswat_goffgratch(dp)  # Goff Gratch formulation, over water

  if (Tconvert is not None):
    eice = esice_goffgratch(dp)  # Goff Gratch formulation, over ice
    ind = num.where(t <= Tconvert)
    e[ind] = eice[ind]
  return e

if __name__ == '__main__':
  print(dp2e.__doc__)
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
  e = dp2e(t,td)
  print(e)
  e = dp2e(t,td,253.15)
  print(e)
