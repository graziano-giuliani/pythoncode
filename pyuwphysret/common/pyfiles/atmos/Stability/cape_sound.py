#!/usr/bin/env python
# cape_sound.py
import numpy as num
from esa import esa

def cape_sound(p,t,r):
  """[cin,cape] = cape_sound(p,t,r)
Calculates convective available potential energy (CAPE) for a sounding. 

Inputs: (arrays)
    p    = atmospheric pressure (mb) at N levels
    t    = dry bulb temperatures (K) at N levels
    r    = water vapor mixing ratio (kg H2O per kg dry air) at N levels

Outputs:
    cape = maximum CAPE for the column (J kg^-1)
    cin  = convective inhibition (J kg^-1)

RLT, 010705"""
  n      = p.size
  tlp    = num.zeros((n,n))
  tlvp   = num.zeros((n,n))
  tvpdif = num.zeros((n,n))
  capep  = num.zeros((n,1))
  nap    = num.zeros((n,1))
  pap    = num.zeros((n,1))

  # thermodynamic constants
  cpd = 1005.7    # specific heat of dry air at const. pressure (J kg-1 K-1)
  cpv = 1870.0    # specific heat of water vapor at const. pressure (J kg-1 K-1)
  cl = 4190.0     # specific heat of liquid water at 0 deg c (J kg-1 K-1)
  cpvmcl = 2320.0 # ?
  rv = 461.5      # gas constant of water vapor = R*/M_H2O (J kg-1 K-1)
  rd = 287.04     # gas constant of dry air = r*/m_dry (J kg-1 K-1)
  eps = rd/rv     # ratio of gas constants (unitless)
  alv0 = 2.501e6  # latent heat of vaporization (J kg-1)
  t0 = 273.15     # C-to-K or K-to-C linear conversion (deg)

  # water vapor pressure ev and saturation vapor pressure es
  tc = t - t0           # temperature in celsius
  ev = r * p/(eps + r)  # vapor pressure
  es = esa(tc)          # saturation vapor pressure

  # begin outer loop, which cycles through parcel origin levels i
  # do calculation only for lowest n/3 levels
  for i in range(0,int(n/3)):
    # define the pseudo-adiabatic entropy sp (conserved variable)
    rs = eps * es[i]/(p[i] - es[i])
    alv = alv0 - cpvmcl * tc[i]
    em = max(ev[i],1e-6)
    sp = (cpd * num.log(t[i]) - rd * num.log(p[i] - ev[i]) +
          alv * r[i] / t[i] - r[i] * rv * num.log(em / es[i]))

    # find lifted condensation pressure plcl
    rh = r[i] / rs  # relative humidity
    rh = max(0.0,min(rh,1.0))
    chi = t[i] / (1669.0 - 122.0 * rh - t[i])
    plcl = p[i] * pow(rh,chi)

    # begin updraft loop
    xsum = 0
    rg0 = r[i]
    tg0 = t[i]
    for j in range(i,n):
      # calculate estimates of the rates of change of the entropies with
      # temperature at constant pressure
      rs = eps * es[j] / (p[j] - es[j])  # saturation mixing ratio
      alv = alv0 - cpvmcl * tc[j]
      slp = (cpd + rs * cl + alv * alv * rs / (rv * t[j] * t[j])) / t[j]
      # calculate lifted parcel temperature below its lcl
      if (p[j] >= plcl):
        tlp[i,j]    = t[i] * pow((p[j] / p[i]),(rd / cpd))
        tlvp[i,j]   = tlp[i,j] * (1.0 + r[i] / eps) / (1.0 + r[i])
        tvpdif[i,j] = tlvp[i,j] - t[j] * (1.0 + r[j] / eps)/(1.0 + r[j])
      else:
        # iteratively calculate lifted parcel temperature and mixing ratios for 
        # pseudo-adiabatic ascent
        tg = t[j]
        rg = rs
        for k in range(0,7):
          cpw = (j>1) * (xsum+cl*0.5*(rg0 + rg)*(num.log(tg) - num.log(tg0)))
          em = rg * p[j] / (eps + rg)
          alv = alv0 - cpvmcl * (tg - 273.15)
          spg = cpd*num.log(tg) - rd*num.log(p[j] - em) + cpw + alv * rg / tg
          tg = tg + (sp - spg) / slp
          enew = esa(tg - 273.15)
          rg = eps * enew / (p[j] - enew)           
        tlp[i,j]    = tg
        tlvp[i,j]   = tg * (1.0 + rg / eps) / (1.0 + rg)
        tvpdif[i,j] = tlvp[i,j] - t[j] * (1.0 + r[j] / eps) / (1.0 + r[j])
        rg0 = rg
        tg0 = tg
        xsum = cpw

    # find positive and negative areas  pa and na and cape (=pa-na) from 
    # pseudo-adiabatic ascent

    # find lifted condensation level and maximum level of positive buoyancy
    #   icb = n  # index of lifted cond level
    inbp = 1 # index of maximum level of positive buoyancy
    for j in range(n-1,i,-1):
#     if (p[j] < plcl):
#       icb = min(icb,j)
      if (tvpdif[i,j] > 0.0):
        inbp = max(j,inbp)
        break
    imax = max(inbp,i)
    tvpdif[i,imax:n-1] = 0  # set to zero for levels above imax
    # do updraft loops
    if (inbp > i):
      for j in range(i+1,inbp+1):
        tvm = 0.5 * (tvpdif[i,j] + tvpdif[i,j-1])
        pm  = 0.5 * (p[j] + p[j-1])
        if (tvm <= 0):
          nap[i] = nap[i] - rd * tvm * (p[j-1] - p[j]) / pm
        else: 
          pap[i] = pap[i] + rd * tvm * (p[j-1] - p[j]) / pm
    capep[i] = pap[i] - nap[i]
  # else cape=0 if no free convection is possible
  cape = max(capep)
  cin  = max(nap)
  return(cape,cin)

if __name__ == '__main__':
  print(cape_sound.__doc__)
  t = num.array(
      ( 24.54, 23.16, 21.67, 20.23, 18.86, 17.49, 16.10, 14.69, 13.22, 11.52,
         9.53,  7.24,  4.80,  2.34,  0.04, -2.29, -4.84, -7.64,-10.66,-13.95,
       -17.54,-21.45,-25.58,-29.90,-34.33,-38.94,-43.78,-48.80,-53.94,-58.79,
       -63.27,-67.32,-70.74,-73.62,-75.74,-77.07,-77.43,-76.63,-75.06,-73.14,
       -71.43 ))
  p = num.array(
      ( 1012.0, 991.3, 969.1, 945.5, 920.4, 893.8, 865.7, 836.1, 805.1, 772.8,
         739.5, 705.2, 670.3, 635.0, 599.7, 564.5, 529.8, 495.7, 462.6, 430.7,
         400.0, 370.8, 343.0, 316.7, 292.0, 266.8, 247.2, 227.0, 208.2, 190.8,
         174.7, 159.9, 146.2, 133.6, 121.9, 111.3, 101.5,  92.6,  84.4,  76.9,
          70.0 ))
  r = num.array(
      ( 17.78, 16.92, 15.93, 14.87, 13.78, 12.70, 11.84, 10.96, 10.15,  9.31,
         8.46,  7.73,  7.05,  6.32,  5.62,  4.91,  4.10,  3.30,  2.67,  2.15,
         1.66,  1.26,  0.95,  0.68,  0.45,  0.28,  0.17,  0.10,  0.06,  0.04,
         0.02,  0.02,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,
         0.02 ))
  t = t + 273.15
  r = r / 1000.0
  [cape,cin] = cape_sound(p,t,r)
  print(cape)
  print(cin)
