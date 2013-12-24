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

#
# Unit test of the above class
#
if ( __name__ == '__main__' ):
  import numpy as np
  a = emissivity('../data/emissivity.nc')
  e = a.get('SfGrd')
  print(e[1])
  print(e[-1])
