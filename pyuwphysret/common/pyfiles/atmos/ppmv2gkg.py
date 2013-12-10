#!/usr/bin/env python
# ppmv2gkg.py
import numpy as num
from gas_mass import gas_mass

def ppmv2gkg(ppmv,gasid):
  """gkg = ppmv2gkg(ppmv,gasid);

function ppmv = ppmv2gkg(ppmv,gasid);

convert volume mixing ratio in ppmv to 
mass mixing ratio in g/kg.

inputs:
    gkg:  mass mixing ratio (g/kg)
  gasid:  HITRAN gas id
  """
  # convert to parts per volume
  ppv = ppmv / 1e6
  # multiply by ratio of masses to get mass 
  # mixing ratio in g/g
  gg = ppv * gas_mass(gasid)/gas_mass(99)
  # multiply by 1000 to get g/kg
  gkg = gg * 1000
  return gkg

if __name__ == '__main__':
  print(ppmv2gkg.__doc__)
  gkg = ppmv2gkg(0.00010718620506252184,23)
  print(gkg)
  try:
    gkg = ppmv2gkg(1e-7,58)
  except:
    print('Ok')
  gkg = ppmv2gkg(6.5773353106547515e-06,'CO2')
  print(gkg)
