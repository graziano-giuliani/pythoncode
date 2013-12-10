#!/usr/bin/env python
# esice_goffgratch.m
import numpy as num

def esice_goffgratch(T):
  """svp = esice_goffgratch(T)

Compute water vapor saturation pressure over ice
using Goff-Gratch formulation.  Adopted from PvD's 
svp_ice.pro.

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
  ewi = 6.1071
  c1 =  9.09718
  c2 = 3.56654
  c3 = 0.876793
  ratio = 273.15/T

  tmp = (( -c1 * ( ratio - 1.0 ) ) -
         (  c2 * num.log10( ratio ) ) +
         (  c3 * ( 1.0 - ( 1.0 / ratio ) ) ) +
         num.log10( ewi ))

  svp = 10.0**tmp
  return svp

if __name__ == '__main__':
  print(esice_goffgratch.__doc__)
  t = num.array(
      ( 24.54, 23.16, 21.67, 20.23, 18.86, 17.49, 16.10, 14.69, 13.22, 11.52,
         9.53,  7.24,  4.80,  2.34,  0.04, -2.29, -4.84, -7.64,-10.66,-13.95,
       -17.54,-21.45,-25.58,-29.90,-34.33,-38.94,-43.78,-48.80,-53.94,-58.79,
       -63.27,-67.32,-70.74,-73.62,-75.74,-77.07,-77.43,-76.63,-75.06,-73.14,
       -71.43 ))
  t = t + 273.15
  e = esice_goffgratch(t)
  print(e)
