#!/usr/bin/env python
# esat.py
import numpy as num

def esat(t):
  """es = esat(t)
Calculates saturation vapor pressure es for given temperature(s) t, using
an adapted form of the Clausius-Clapeyron equation. 

Inputs:
    t    = dry bulb temperature(s) (deg C or K) (array or numeric)

Outputs:
    esat = saturation vapor pressure (mb) for the given temperature

RLT, 010607"""
  e1 = 1013.250
  tk = 273.16
  try:
    if (max(t) < 105.0):
      t0 = tk
    else:
      t0 = 0.0
  except:
    if (t < 105.0):
      t0 = tk
    else:
      t0 = 0.0
  # Required temeprature in K
  esat = e1 * pow(10.0,10.79586*(1.0-tk/(t+t0))-5.02808*num.log10((t+t0)/tk) +
                  1.50474*1.0e-4*(1.0-pow(10.0,-8.29692*((t+t0)/tk-1.0))) +
                  0.42873*1.0e-3*(pow(10.0,4.76955*(1.0-tk/(t+t0)))-1.0) -
                  2.2195983)
  return esat

if __name__ == '__main__':
  print(esat.__doc__)
  t = num.array(
      ( 24.54, 23.16, 21.67, 20.23, 18.86, 17.49, 16.10, 14.69, 13.22, 11.52,
         9.53,  7.24,  4.80,  2.34,  0.04, -2.29, -4.84, -7.64,-10.66,-13.95,
       -17.54,-21.45,-25.58,-29.90,-34.33,-38.94,-43.78,-48.80,-53.94,-58.79,
       -63.27,-67.32,-70.74,-73.62,-75.74,-77.07,-77.43,-76.63,-75.06,-73.14,
       -71.43 )) + 273.15
  es = esat(t)
  print(es)
