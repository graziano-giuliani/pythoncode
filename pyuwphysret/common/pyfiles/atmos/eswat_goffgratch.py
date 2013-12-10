#!/usr/bin/env python
# eswat_goffgratch.py
import numpy as num

def eswat_goffgratch(T):
  """svp = eswat_goffgratch(T)

Compute water vapor saturation pressure over water
using Goff-Gratch formulation.  Adopted from PvD's 
svp_water.pro.

Inputs: 
   T       temperature [Kelvin]

Output:
   svp    saturation pressure [mbar]

Notes: svp returned for all values of input T,
   but results not valid for T >= 370 K and 
   T <= 160 K.

Reference: Goff-Gratch formulation from sixth revised 
      edition of Smithsonian Meteorology Tables.

DCT 8/22/00
  """
  t_sat = 373.16
  t_ratio = t_sat/T
  rt_ratio = 1.0/t_ratio
  sl_pressure = 1013.246

  c1 = 7.90298
  c2 = 5.02808
  c3 = 1.3816e-7
  c4 = 11.344
  c5 = 8.1328e-3
  c6 = 3.49149

  tmp = (( -1.0 * c1 * ( t_ratio - 1.0 ) ) +
         ( c2 * num.log10( t_ratio ) ) -
         ( c3 * ( 10.0**( c4 * ( 1.0 - rt_ratio ) ) - 1.0 ) ) +
         ( c5 * ( 10.0**( -1.0 * c6 * ( t_ratio - 1.0 ) ) - 1.0 ) ) +
         num.log10( sl_pressure ))
  svp = 10.0**tmp;
  return svp

if __name__ == '__main__':
  print(eswat_goffgratch.__doc__)
  t = num.array(
      ( 24.54, 23.16, 21.67, 20.23, 18.86, 17.49, 16.10, 14.69, 13.22, 11.52,
         9.53,  7.24,  4.80,  2.34,  0.04, -2.29, -4.84, -7.64,-10.66,-13.95,
       -17.54,-21.45,-25.58,-29.90,-34.33,-38.94,-43.78,-48.80,-53.94,-58.79,
       -63.27,-67.32,-70.74,-73.62,-75.74,-77.07,-77.43,-76.63,-75.06,-73.14,
       -71.43 ))
  t = t + 273.15
  e = eswat_goffgratch(t)
  print(e)
