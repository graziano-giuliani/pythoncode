#!/usr/bin/env python
# mr2e.py
import numpy as num
from esice_goffgratch import esice_goffgratch
from eswat_goffgratch import eswat_goffgratch

def mr2e(p,mr):
  """e = mr2e(p,mr)

compute H2O partial pressure (e,mbar) given
pressure (p,mbar) and H2O mass mixing ratio (mr,g/kg)

DCT 3/6/00
  """
  # ratio of water mass to dry air mass
  eps = 0.621970585
  e = p*mr/(1000*eps + mr)
  return e

if __name__ == '__main__':
  print(mr2e.__doc__)
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
  e = mr2e(p,r)
  print(e)
