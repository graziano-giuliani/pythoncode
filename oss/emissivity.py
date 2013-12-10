#!/usr/bin/env python

class emissivity:
  """Atmospheric emissivity dataset"""
  def __init__(self,datafile):
    from netCDF4 import Dataset
    self.df = Dataset(datafile)
  def get(self,f,part=None):
    import numpy as np
    if ( part is None ):
      return np.array(self.df.variables[f])
    else:
      return np.array(self.df.variables[f][part])
  def __del__(self):
    self.df.close()

if ( __name__ == '__main__' ):
  import numpy as np
  a = emissivity('../data/emissivity.nc')
  e = a.get('SfGrd')
  print(e[1])
  print(e[-1])
