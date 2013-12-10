#!/usr/bin/env python

class hitran:
  """Loads hitran pre-computed tables in memory"""
  def __init__(self,datafile):
    from netCDF4 import Dataset
    self.df = Dataset(datafile)
    self.desc = self.df.info
    self.v1 = self.df.v1
    self.v2 = self.df.v2
  def get(self,name,part=None):
    import numpy as np
    if ( part is None ):
      return np.array(self.df.variables[name][:])
    else:
      return np.array(self.df.variables[name][part])
  def nchan(self):
    return len(self.df.dimensions['nchan'])
  def nlayod(self):
    return len(self.df.dimensions['nlayod'])
  def ntmpod(self):
    return len(self.df.dimensions['ntmpod'])
  def nf(self):
    return len(self.df.dimensions['nf'])
  def nmol(self):
    return len(self.df.dimensions['nmol'])
  def mxmols(self):
    return len(self.df.dimensions['mxmols'])
  def nmoltab(self):
    return len(self.df.dimensions['nmoltab'])
  def __del__(self):
    self.df.close( )

if ( __name__ == '__main__' ):
  import numpy as np
  a = hitran('../data/leo.iasi.0.05.nc')
  print(a.desc)
  wvptab = a.get('wvptab')
  imols = a.get('imols')
  inhmap = a.get('ichmap_ir',[range(0,20),range(0,10)])
  print(np.shape(inhmap))
  inhmap = a.get('ichmap_ir')
  print(np.shape(inhmap))
  print(a.nchan())
