#!/usr/bin/env python
# satvap.py
import numpy as num
from esice_goffgratch import esice_goffgratch
from eswat_goffgratch import eswat_goffgratch

def satvap(T,Tconvert=None):
  """esat = satvap(T,Tconvert)

compute saturation vapor pressure [mbar] given temperature, T [K].
If Tconvert is input, the calculation uses the saturation vapor 
pressure over ice (opposed to over water) for temperatures less than 
Tconvert [K].

DCT, updated 3/5/00
  """
  #saturation pressure over water
  esat = eswat_goffgratch(T) # Goff Gratch formulation, over water

  # saturation pressure over ice if needed
  if (Tconvert is not None):
    try:
      ind = num.where(T <= Tconvert)
      esat[ind] = esice_goffgratch(T[ind]); # Goff Gratch formulation, over ice
    except:
      if (T <= Tconvert):
        esat = num.array((esice_goffgratch(T)))
  return esat

if __name__ == '__main__':
  print(satvap.__doc__)
  t = num.array(
      ( 24.54, 23.16, 21.67, 20.23, 18.86, 17.49, 16.10, 14.69, 13.22, 11.52,
         9.53,  7.24,  4.80,  2.34,  0.04, -2.29, -4.84, -7.64,-10.66,-13.95,
       -17.54,-21.45,-25.58,-29.90,-34.33,-38.94,-43.78,-48.80,-53.94,-58.79,
       -63.27,-67.32,-70.74,-73.62,-75.74,-77.07,-77.43,-76.63,-75.06,-73.14,
       -71.43 ))
  t = t + 273.15
  e = satvap(t,253.15)
  print(e)
