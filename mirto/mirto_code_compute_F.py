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
    self.cx = control
    self.xdim = control.xdim
    self.outdata = {}
    self.oss = control.oss
    self.wnF = self.oss.cwvn
    self.F = np.zeros(len(self.oss.cwvn))
    self.Md = 28.966  # Molecular mass of dry air
    self.Mw = 18.016  # Molecular Mass of water
    self.Mc = 44.01   # Molecular Mass of CO2
    self.Mo = 48      # Molecular Mass of Ozone
    self.tds = 0
    self.tde = self.tds+self.xdim[0]
    self.wvs = self.tde
    self.wve = self.wvs+self.xdim[1]
    self.cos = self.wve
    self.coe = self.cos+self.xdim[2]
    self.ozs = self.coe
    self.oze = self.coe+self.xdim[3]
    self.sks = self.oze
    self.ske = self.sks+self.xdim[4]
    self.ems = self.ske
    self.eme = self.ems+self.xdim[5]

  def compute(self,x):
    from mirto_code_configuration import surface_emissivity
    # Define Input Parameters
    # Construct prof structure from input state vector (x, p)
    # Note: The state vector unit for concentration is the natural log
    # of vmr in ppv. 
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

    vmr = 1.E6*np.exp(x[self.wvs:self.wve])
    co2_ppmv = x[self.cos:self.coe]
    oz_vmr = x[self.ozs:self.oze]

    indata = {}
    emiss = surface_emissivity(self.cx,x)
    indata['sfgrd'] = emiss.wnSurfEmiss
    indata['emrf'] = emiss.SurfEmiss_values
    indata['tskin'] = x[self.sks:self.ske]
    indata['psf'] = self.cx.surfacePressure_mb
    indata['temp'] = np.flipud(x[self.tds:self.tde])
    indata['h2o'] = np.flipud(1.E-6*vmr)
    indata['co2'] = np.flipud(1.E-6*(self.Mc/self.Md)*co2_ppmv)
    indata['o3'] = np.flipud(np.exp(oz_vmr))
    indata['pressure'] = np.flipud(self.cx.pressure_grid)
    indata['pobs'] = self.cx.observationPressure_mb
    indata['obsang'] = self.cx.FOVangle
    indata['sunang'] = self.cx.SunAngle
    self.oss.compute(indata,self.outdata)
    self.F = self.outdata['y']

  def estimate_K(self,x):
    """
Compute the jacobian K using the selected model
    """
    xdim = self.xdim
    Jvar = self.cx.Jvar
    SEflag = np.count_nonzero(np.where(Jvar==-2))
    SKTflag = np.count_nonzero(np.where(Jvar==-1))
    Tflag = np.count_nonzero(np.where(Jvar==0))
    WVflag = np.count_nonzero(np.where(Jvar==1))
    CO2flag = np.count_nonzero(np.where(Jvar==2))
    O3flag = np.count_nonzero(np.where(Jvar==3))
    
    # Set number of channels
    krow = len(self.cx.wnR)

    # Inizialize the K matrix with nans
    jac = np.zeros((krow,sum(xdim[0:5])))
    jac.fill(np.NAN)
    wvn = self.cx.wnR

    if ( Tflag > 0 ):
      jcb = np.transpose(self.outdata['xkt'][self.tds:self.tde,:])
      jac[:,self.tds:self.tde] = np.fliplr(jcb)
    if ( WVflag > 0 ):
      #  convert from log(vmr in ppv) to vmr in ppmv
      vmr = 1.e6*np.exp(x[self.wvs:self.wve])
      # w = 1.e-6*(self.Mw/self.Md)*vmr 
      w = 1.e-6*vmr
      w_mat = np.tile(w,(krow,1))
      jcb = np.transpose(self.outdata['xkt'][self.wvs:self.wve,:])
      # Jacobians in log(q)
      jac[:,self.wvs:self.wve] = np.fliplr(jcb)*w_mat
    if ( CO2flag > 0 ):
      # convert from log(vmr in ppv) to vmr in ppmv
      # vmr = x[self.cos:self.coe]
      # w = 1.e-6*(self.Mc/self.Md)*vmr
      # w_mat = np.tile(w,(krow,1))
      jcb = np.transpose(self.outdata['xkt'][self.cos:self.coe,:])
      # jac[:,self.cos:self.coe] = np.fliplr(jcb)*w_mat
      # Jacobians in ppmv
      jac[:,self.cos:self.coe] = np.fliplr(jcb)*(self.Mc/self.Md)*1.e-6
    if ( O3flag > 0 ):
      # convert from log(vmr in ppv) to vmr in ppmv
      vmr = np.exp(x[self.ozs:self.oze])
      # w = 1.e-6*(self.Mo/self.Md)*vmr
      w = vmr
      w_mat = np.tile(w,(krow,1))
      jcb = np.transpose(self.outdata['xkt'][self.ozs:self.oze,:])
      # jacobians in log(q)
      jac[:,self.ozs:self.oze] = np.fliplr(jcb)*w_mat
    if ( SKTflag > 0 ):
      jcb = np.transpose(self.outdata['xkt'][self.sks:self.ske,:])
      jac[:,self.sks:self.ske] = jcb
    if ( SEflag > 0 ):
      # margin used in calculation about selwn
      calcmargin = 70
      # compute surface emissivity
      emiss = surface_emissivity(self.cx,x)
      interp_EmissVal = np.interp(wvn,
            emis.wnSurfEmiss,emis.SurfEmiss_values,0.999,0.999)
      jcb = np.transpose(self.outdata['paxkemrf'][1,:])
      self.cx.Kse = jcb
      w = interp_EmissVal*(1.0-interp_EmissVal)
      self.cx.Kse_logit = jcb*w
      # Note: The jacobian we want is the lblrtm surface jacobian times each
      # SurfEmissModelFunctions interpolated to the calculation scale.
      interp_SurfEmissModelFunctions = np.interp(wvn,
            emis.wnSurfEmissModelFunctions,emis.SurfEmissModelFunctions,0.0,0.0)
      for i in range(self.ems,self.eme):
        jac[:,i] = interp_SurfEmissModelFunctions[:,i]*jcb*w
    self.K = jac
    self.wnK = wvn

#
# Unit test of the above
#
if ( __name__ == '__main__' ):
  from netCDF4 import Dataset
  sys.path.append('/home/graziano/Software/pythoncode/oss')
  from oss4SHIS import oss4SHIS
  #
  # OSS init input
  #
  solar = 'solar_irradiances.nc'
  precomputed = 'leo.cris.0.05.nc'
  datapath = '/home/graziano/Software/pythoncode/data'
  oss = oss4SHIS(os.path.join(datapath,solar),
                 os.path.join(datapath,precomputed))
  from mirto_code_configuration import control , apriori
  cx = control(datapath,oss)
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
  rad.estimate_K(x)
  print(rad.K)
