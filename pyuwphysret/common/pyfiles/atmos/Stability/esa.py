#!/usr/bin/env python
# esa.py
import numpy as num

def esa(t):
  """es = esa(t)
Calculates saturation vapor pressure es for given temperature(s) t, using
an adapted form of the Clausius-Clapeyron equation. 
 
Inputs:
    t  = dry bulb temperature(s) (deg C o K) (array or numeric)

Outputs:
    es = saturation vapor pressure (mb) for the given temperature

RLT, 010530
  """
  try:
    if (max(t) > 105.0):
      t0 = 273.15
    else:
      t0 = 0.0
  except:
    if (t > 105.0):
      t0 = 273.15
    else:
      t0 = 0.0
  es = 6.112*num.exp(17.67*(t-t0)/(243.5+(t-t0)))
  return es

if __name__ == '__main__':
  print(esa.__doc__)
  t = num.linspace(300.0,200.0,num=20)
  es = esa(t)
  print(t)
  print(es)
  tc = 300.0-273.15
  es = esa(tc)
  print(tc,es)
