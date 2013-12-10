#!/usr/bin/env python

class asolar:
  """Solar irradiance dataset"""
  def __init__(self,datafile):
    from netCDF4 import Dataset
    import numpy as np
    df = Dataset(datafile)
    self.frq = np.array(df.variables['FREQ'][:],np.double)
    self.irr = np.array(df.variables['IRRADIANCE'][:],np.double)
    df.close()
  def get(self,f):
    import numpy as np
    i = np.interp(f,self.frq,self.irr)
    return i

if ( __name__ == '__main__' ):
  import numpy as np
  a = asolar('../data/solar_irradiances.nc')
  for j in range(1,1000):
    f =  np.linspace(1550+j, 2550+j+1, 1550)
    i = a.get(f)
    print(j)
    print(f[0],f[1549])
    print(i[0],i[1549])
