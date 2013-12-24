#!/usr/bin/env python
# mirto_code_compute_F.py
import numpy as np
import os
import sys

class radiance:
  """radiance class

This function is part of the Python implementation of the MTG-IRS Real Time
Operational code

Compute the radiance vector F using the selected forward model

control is a data structure containing control variables
x is a state vector containing the variables for which the forward
model should be evaluated

Note: The state vector unit for concentration is the natural log or vmr in ppv. 
xdim is a vector of dimensions for each of parameters of the solution
 (T, wv, ozone, surface)

Calling sequence:

Developed by: P. Anontelli, G. Giuliani, T. Cherubini

Wed 11 dec 2013, 15.59.41, CET
  """
  def __init__(self,control):
    #  fm.wnF = calculation wavenumber vector
    #  fm.F   = calculation radiance vector
    sys.path.append(control.osspath)
    from oss4SHIS import oss4SHIS
    self.cx = control
    self.outdata = {}
    self.oss = oss4SHIS(os.path.join(control.datapath,control.solar),
                        os.path.join(control.datapath,control.precomputed))
    self.wnF = self.oss.cwvn
    self.F = np.zeros(len(self.oss.cwvn))
  def compute(self,x):
    from mirto_code_configuration import surface_emissivity
    # Define Input Parameters
    # Construct prof structure from input state vector (x, p)
    # Note: The state vector unit for concentration is the natural log
    # of vmr in ppv. 
    Md = 28.966  # Molecular mass of dry air
    Mw = 18.016  # Molecular Mass of water
    Mc = 44.01   # Molecular Mass of CO2
    Mo = 48      # Molecular Mass of Ozone
    #
    # Define Input to radiance calculator
    #   Atmospheric profile data: 
    #      pressure (Nx1) level pressure (mbar)
    #      tdry     (Nx1) level temperature (K)
    #      w        (Nx1) level water vapor mass mixing ratio (g/kg)
    #      vmr      (Nx1) level water vapor volume mixing ratio (ppmv)
    #      alt      (Nx1) level altitudes (km)   
    #      oz       (Nx1) level ozone (ppmv)
    #
    # NOTE: when using pressure boundaries the surface elevation is required
    # in the lowest altitude level for use in the hypsometric equation.
    # All other altitudes are ignored in the profile level input. 
    # 
    # convert from log(vmr in ppv) to vmr in ppmv
    xdim = self.cx.xdim

    tds = 0
    tde = tds+xdim[0]
    wvs = tde
    wve = wvs+xdim[1]
    cos = wve
    coe = cos+xdim[2]
    ozs = coe
    oze = coe+xdim[3]
    sk = oze

    vmr = 1.E6*np.exp(x[wvs:wve])
    co2_ppmv = x[cos:coe]
    oz_vmr = x[ozs:oze]

    indata = {}
    emiss = surface_emissivity(self.cx,x)
    indata['sfgrd'] = emiss.wnSurfEmiss
    indata['emrf'] = emiss.SurfEmiss_values
    indata['tskin'] = x[sk]
    indata['psf'] = self.cx.surfacePressure_mb
    indata['temp'] = np.flipud(x[tds:tde])
    indata['h2o'] = np.flipud(1.E-6*vmr)
    indata['co2'] = np.flipud(1.E-6*(Mc/Md)*co2_ppmv)
    indata['o3'] = np.flipud(np.exp(oz_vmr))
    indata['pressure'] = np.flipud(self.cx.pressure_grid)
    indata['pobs'] = self.cx.observationPressure_mb
    indata['obsang'] = self.cx.FOVangle
    indata['sunang'] = self.cx.SunAngle
    print (indata)
    self.oss.compute(indata,self.outdata)
    self.F = self.outdata['y']

#
# Unit test of the above
#
if ( __name__ == '__main__' ):
  from netCDF4 import Dataset
  from mirto_code_configuration import control , apriori
  cx = control('/home/graziano/Software/pythoncode/data',
               '/home/graziano/Software/pythoncode/oss')
  ap = apriori(cx)
  rad = radiance(cx)
  x = ap.x0
  rad.compute(x)
  rootgrp = Dataset('radiance.nc', 'w', format='NETCDF4')
  nchan = rootgrp.createDimension('nchan', len(rad.wnF))
  ncwn = rootgrp.createVariable('wn','f4',('nchan',))
  ncra = rootgrp.createVariable('radiance','f4',('nchan',))
  ncwn[:] = rad.wnF
  ncra[:] = rad.F
  rootgrp.close()
