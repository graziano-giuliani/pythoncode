#!/usr/bin/env python
# gkg2ppmv.py
import numpy as num
from gas_mass import gas_mass

def gkg2ppmv(gkg,gasid):
  """ppmv = gkg2ppmv(gkg,gasid)

convert mass mixing ratio in g/kg to 
volume mixing ratio in ppmv.

inputs:
    gkg:  mass mixing ratio (g/kg)
  gasid:  HITRAN gas id
  """
  #divide by 1000 to get g/g
  gg = gkg / 1000

  # divide by ratio of masses to get volume
  # mixing ratio in g/g
  ppv = gg / (gas_mass(gasid)/gas_mass(99))

  # convert to parts per million volume
  ppmv = ppv * 1e6
  return ppmv

if __name__ == '__main__':
  print(gkg2ppmv.__doc__)
  ppmv = gkg2ppmv(1e-7,23)
  print(ppmv)
  try:
    mm = gkg2ppmv(1e-7,58)
  except:
    print('Ok')
  ppmv = gkg2ppmv(1e-8,'CO2')
  print(ppmv)
