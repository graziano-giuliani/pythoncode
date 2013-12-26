#!/usr/bin/env python
# mirto_code_configuration.py
"""
This module is a temporary for a test configuration of the mirto code
"""
import numpy as np
from netCDF4 import Dataset
import os

class mirto_config:
  """mirto_config

This function is part of the Python implementation of the MTG-IRS Real Time
Operational code

Defines input variable for test case for development of python version of
UWPHYSRET. Test case selected for IASI radiance spectrum.

Calling sequence:

Developed by: P. Anontelli, G. Giuliani, T. Cherubini

Wed 11 dec 2013, 15.59.41, CET
  """
  def __init__(self,datapath,oss):
    #
    # Surface Elevation
    #
    # self.surfaceElevation_km = 0.0145
    self.datapath = datapath
    self.oss = oss
    self.surfacePressure_mb = 1013.0
    self.pressure_surf = self.surfacePressure_mb
    #
    # Satellite/FOV properties
    #
    self.observationAltitude_km = 100.0 # not used
    self.observationPressure_mb = 0.005
    self.SunAngle = 90.0
    #
    # Define other retrieval parameters
    #
    self.gamma = 3.0  # Marquardt-Levemberg parameter
    # Set a config flag(0=Use existing K or 1=Update K) for each desired
    # iteration (including the first)
    # Note: Updating the jacobian matrix K is very slow and may be unnecessary
    # Enter config flags for first retrieval,      e.g. [1 0 0 0 0 0]
    self.EstimateK = np.array([1,1,1,1,1,1,1,1,1,1,1],dtype=int)
    # Enter Environmental Parameter to Retrieve (e.g. [-2 -1 0 1 2 3])
    self.Jvar = np.array([-2,-1,0,1,3],dtype=int)
    self.Iteration_limit = sum(self.EstimateK)
    # Get input parameters
    # Array of observations
    df = Dataset(os.path.join(self.datapath,'fov.nc'))
    self.radiance = np.array(df.variables['selradiances'][:])
    # Field of View Angle
    self.FOVangle = np.array(df.variables['FOVangle'][:])
    # Indeces for selected channels to be used with Forward Model
    # Python indexes start from 0
    self.indx = np.array(df.variables['indxselchannel'][:],dtype=int)-1
    # FOV latitude
    self.fov_latitude = np.array(df.variables['Latitude'][:])
    # whole observation wavenumber grid
    self.wnR = np.array(df.variables['Wavenumber'][:])
    # whole observation radiance spectrum
    self.R = np.array(df.variables['Radiance'][:])
    df.close( )
    #
    df = Dataset(os.path.join(self.datapath,'fg.nc'))
    self.p = np.array(df.variables['p'][:])
    self.xdim = np.array(df.variables['xdim'][:],dtype=int)
    df.close( )
    # Get emissivity Principal Componets (model functions)
    df = Dataset(os.path.join(self.datapath,'emissivitymodel.nc'))
    self.SurfEmissModelFunctions = np.array(df.variables['ModelFunctions'][:])
    self.SurfEmissModelFunctions_bias = np.array(
            df.variables['ModelFunctionsBias'][:])
    self.wnSurfEmissModelFunctions = np.array(
            df.variables['WnModelFunctions'][:])
    df.close()
    # Get indices for selected variables within whole state vector
    df = Dataset(os.path.join(self.datapath,'apriori.nc'))
    self.state_var_indx = np.array(df.variables['varindx'][:],dtype=int)-1
    df.close()
    self.pressure_grid = self.p[0:self.xdim[0]]

class mirto_obsErr:
  def __init__(self,config):
    df = Dataset(os.path.join(config.datapath,'obserr.nc'))
    # Get obervation error covariance
    self.Se = np.array(df.variables['obserrselchannels'][:])
    self.SeInv = np.array(df.variables['invobserrselchannels'][:])
    df.close()

class mirto_apriori:
  def __init__(self,config):
    df = Dataset(os.path.join(config.datapath,'fg.nc'))
    # Get first guess
    self.x0 = np.array(df.variables['x0'][:])
    self.xa = np.array(df.variables['xa'][:])
    df.close()
    df = Dataset(os.path.join(config.datapath,'apriori.nc'))
    # Get a-priori covariance
    self.Sa = np.array(df.variables['Sa'][:])
    self.SaInv = np.array(df.variables['SaInv'][:])
    df.close()

class surface_emissivity:
  def __init__(self,config,x):
    """
Compute a surface emissivity spectrum from a set of model coefficients and
model functions

   out1 = computed emissivity values
   out2 = wavenumber scale of computed emissivity values

   config is a data structure containing config variables
   SurfEmissCoefficients are the surface emissivity coefficients
       (vector)
   SurfEmissModelFunctions are the surface emissivity model functions
       (matrix with functions in rows)
   wnSurfEmissModelFunctions is the wavenumber scale of the surface
       emissivity model funtions (vector)
   SurfEmissModelFunctions_bias is the offset (bias) for all the model
       functions
   wnSurfEmissModelFunctions_bias is the wavenumber scale for the bias
       spectrum.
   OutputStartWn is the desired output starting wavenumber
   OutputDeltaWn is the desired output wavenumber interval
   OutputNumberPoints is the desired length of the output emissivity
       spectrum

Note:
 The surface emissivity model functions are defined on their own
 wavenumber scale.
 This routine handles any interpolation or extrapolation to the requested
 output range.
 This routine does not perform any truncation of output values, this is
 left to the calling program.

Calling sequence:

Developed by: P. Anontelli, G. Giuliani, T. Cherubini

Wed 11 dec 2013, 15.59.41, CET
    """
    x_dim = config.xdim
    # load surface emissivity coefficients from solution vector
    SurfEmissCoefficients = (
          x[x_dim[0]+x_dim[1]+x_dim[2]+x_dim[3]+x_dim[4]:
            x_dim[0]+x_dim[1]+x_dim[2]+x_dim[3]+x_dim[4]+x_dim[5]] )
    SurfEmissModel_values = ( np.dot(config.SurfEmissModelFunctions.T,
      SurfEmissCoefficients) + config.SurfEmissModelFunctions_bias)
    wnSurfEmiss = config.wnSurfEmissModelFunctions
    self.SurfEmiss_values = SurfEmissModel_values
    self.wnSurfEmiss = wnSurfEmiss
    self.SurfEmiss_values = ((    np.exp(self.SurfEmiss_values))/
                             (1.0+np.exp(self.SurfEmiss_values)))
    self.SurfEmissModelFunctions = config.SurfEmissModelFunctions
    self.wnSurfEmissModelFunctions = config.wnSurfEmissModelFunctions
    self.SurfEmissModelFunctions_bias = config.SurfEmissModelFunctions_bias
    self.wnSurfEmissModelFunctions_bias = config.wnSurfEmissModelFunctions

#
# Unit test of the above
#
if ( __name__ == '__main__' ):
  # Do not need oss to test configuration
  datapath = '/home/graziano/Software/pythoncode/data'
  oss = None
  cx = mirto_config(datapath,oss)
  print(cx.pressure_grid)
  oe = mirto_obsErr(cx)
  ap = mirto_apriori(cx)
  x = ap.x0
  emiss = surface_emissivity(cx,x)
  rootgrp = Dataset('emiss.nc', 'w', format='NETCDF4')
  nspe = rootgrp.createDimension('nspectral', len(emiss.wnSurfEmiss))
  ncSfGrd = rootgrp.createVariable('SfGrd','f4',('nspectral',))
  ncEmRf = rootgrp.createVariable('EmRf','f4',('nspectral',))
  ncSfGrd[:] = emiss.wnSurfEmiss
  ncEmRf[:] = emiss.SurfEmiss_values
  rootgrp.close()
