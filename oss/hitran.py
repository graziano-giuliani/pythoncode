#!/usr/bin/env python
#
# Copyright (c) 2013 Paolo Antonelli, Tiziana Cherubini, Graziano Giuliani
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
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

#
# Unit test of the above class
#
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
