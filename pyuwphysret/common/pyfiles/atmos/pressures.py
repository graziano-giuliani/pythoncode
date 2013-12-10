#!/usr/bin/env python
# pressures.m
import numpy as num

def pressures(wflag=1,pmin=0.01,pmax=1100.0,deltalnp=0.06):
  """pres = pressures(wflag,pmin,pmax,deltalnp)

compute level and layer pressures

inputs:
    wflag:  flag deterimining which layering to use
  (defaults to AIRS levels/layers if nargin == 0)

   wflag = 1: use AIRS/kCARTA values
         = 2: pmax to pmin with delta(log(p_level)) = deltalnp
              pmax, pmin, deltalnp default to 1100, 0.01, and 0.06
         = 3: use AERI 100 layers

outputs:  (all in mbar)
    p_level:       level pressure values (nlevels x 1)
    p_level_lower: lower level (wrt altitude) pressure (nlevels-1 x 1)
    p_level_upper: upper level (wrt altitude) pressure (nlevels-1 x 1)
    p_layer:       mean layer pressure (nlevels-1=nlayers x 1)


        ---------------------- p_level_upper(i+1)
                 p_layer(i+1)
        ---------------------- p_level_upper(i) = p_level_lower(i+1)
                 p_layer(i)
        ---------------------- p_level_lower(i)


       ////////////////Earth///////////////////
  """
  if (wflag == 1):
    # Compute AIRS/kCARTA pressures (mbar).  taken from kLAYERS documentation.
    i = num.reshape(num.arange(1,102),(-1,1))
    # constants determined st p_level(1)=1100,p_level(38)=300,p_level(101)=5e-3
    A = -1.5508E-4
    B = -5.5937E-2
    C = 7.4516

    # 101 level pressures
    p_level = num.exp((7.0/2.0)*num.log(A*i**2 + B*i + C))
    # lower and upper (in altitude) levels for each layer

    p_level_lower = p_level[0:100]
    p_level_upper = p_level[1:101]

    # mean layer pressures
    p_layer = ((p_level_lower-p_level_upper) / 
                 num.log(p_level_lower/p_level_upper))
  elif (wflag == 2):
    # levels from pmax to pmin with delta(log(p)) = deltalnp
    lnp = num.zeros((1))
    p_level = num.zeros((1))
    lnpp = num.zeros((1))
    pp_level = num.zeros((1))
    lnp[0] = num.log(pmax)
    p_level[0] = num.exp(lnp[0])
    pp_level[0] = p_level[0]
    while (pp_level[0] >= pmin):
      lnpp[0] = lnp[-1]-deltalnp
      pp_level[0] = num.exp(lnpp[0])
      lnp = num.append(lnp,lnpp,axis=0)
      p_level = num.append(p_level,pp_level,axis=0)
    p_level = num.reshape(p_level,(-1,1))
  elif (wflag == 3):
    # add in 101 AERI retrieval levels here
    pass
  nlevels = p_level.size
  p_level_lower = p_level[0:nlevels-1]
  p_level_upper = p_level[1:nlevels]
  p_layer = (p_level_lower-p_level_upper)/num.log(p_level_lower/p_level_upper)
  pmin = p_level[nlevels-1][0]
  pmax = p_level[0][0]

  # put output variables into a structure
  pres = { }
  pres['p_level'] = p_level
  pres['p_level_lower'] = p_level_lower
  pres['p_level_upper'] = p_level_upper
  pres['p_layer'] = p_layer
  pres['pmin'] = pmin
  pres['pmax'] = pmax
  return pres

if __name__ == '__main__':
  import pprint
  print(pressures.__doc__)
  pp = pprint.PrettyPrinter(indent=4)
  pres = pressures(1)
  pp.pprint(pres)
  pres = pressures(2)
  pp.pprint(pres)
