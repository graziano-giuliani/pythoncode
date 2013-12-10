#!/usr/bin/env python
# rad2bt.py
import numpy as num

def rad2bt(freq,radiance):
  """bt = rad2bt(freq,radiance)

compute brightness temperature given wavenumbers and radiances.

Inputs:
  freq      wavenumbers (Nx1) in 1/cm 
  radiance  radiances (Nx1) in mW/(m^2 sr. 1/cm)
Outputs:
  bt        computed brightnes temperature (Nx1) in Kelvin

see also RADTOT.M, BT2RAD.M, TTORAD.M

DCT 11/11-99

DCT 11-28-05  Output bts are set to NaN if input radiances are < 0
  """
  # fundamental constants:
  #  (Cohen, E.R. and B.N. Taylor, The 1986 CODATA recommended values
  #  of the fundamental physical constants, Journal of Research of
  #  the National Bureau of Standard, 92(2), March-April 1987.)

  h = 6.6260755E-34  # Planck constant in Js
  c = 2.99792458E8   # photon speed in m/s
  k = 1.380658E-23   # Boltzmann constant in J/K

  c1 = 2*h*c*c*1e8
  c2 = h*c/k*1e2

  bt = c2*freq/(num.log(1.0+(c1*freq*freq*freq)/(radiance/1e3)))

  ind = num.where(radiance < 0)
  if ( ind[0].size > 0 ):
    print('WARNING: negative input RAD to RAD2BT, Setting BTs to NaN')
    bt[ind] = num.nan
  return bt

if __name__ == '__main__':
  print(rad2bt.__doc__)
  freq = 400.0
  radiance = 101.0
  bt = rad2bt(freq,radiance)
  print(bt)
