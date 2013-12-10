#!/usr/bin/env python
# stability.py
import numpy as num
from mixr import mixr
from pressure_at_lcl import pressure_at_lcl
from parcel_thetae import parcel_thetae
from temp_along_satadiabat import temp_along_satadiabat
from cape_sound import cape_sound
from tpwv import tpwv

def stability(P,T,Td):
  """StabilityIndices = stability(P,T,Td)
Calculates stability indicies for a sounding. 

Inputs: (all vectors of length N)
    P = atmospheric pressure (mb) at N levels
    T = dry bulb temperatures (K) at N levels
   Td = dewpoint temperatures (K) at N levels

Outputs:
  StabilityIndices = struct array containing the following stability indices:
    LiftedIndex   = Lifted index (degrees C)
    TotalTotals   = Total totals (degrees C)
    CrossTotals   = Cross totals (degrees C)
    VertTotals    = Vertical totals (degrees C)
    KIndex        = K Index (degrees C)
    LCL           = lifted condensation level (mb)
    LFC           = level of free convection (mb)
    Showalter     = Showalter Index (degrees C)
    EqLevel       = equilibrium level (mb)
    CAPE          = maximum CAPE for the column (J kg^-1)
    CIN           = convective inhibition (J kg^-1)
    Parcel_ThetaE = parcel theta-e (K), for a parcel with weighted average 
                    T and Td of lowest 100 mb raised from 50 mb to 500 mb
    TPWV          = total precipitable water vapor (cm)

RLT, 010713
  """
  StabilityIndices = {}
  # Ensures that soundings are properly oriented
  sizeP = num.shape(P)
  try:
    if (sizeP[0] < sizeP[1]):
      P  = num.transpose(P,(-1,1))
      T  = num.transpose(T,(-1,1))
      Td = num.transpose(Td,(-1,1))
  except:
    pass

  if (P[0] < 800.0):
    P  = flipud(P)
    T  = flipud(T)
    Td = flipud(Td)

  # N = number of levels in the sounding
  N = P.size

  # Kelvin temperature check
  if (min(T) < 105):
    T = T + 273.15
  if (min(Td) < 105):
    Td = Td + 273.15

  # Computes mixing ratio in g H20/kg air
  R = mixr(Td,P)/1000.0;

  # Finds levels closest to 850, 700, 500 and 300 mb.
  N_850 = num.where(abs(P-850) == min(abs(P-850)))
  N_700 = num.where(abs(P-700) == min(abs(P-700)))
  N_500 = num.where(abs(P-500) == min(abs(P-500)))
  N_300 = num.where(abs(P-300) == min(abs(P-300)))

  #
  # *** Lifted Index (degrees C) ***
  # LiftedIndex = T_500 - T_parcel
  # T_parcel = parcel with weighted average T and Td of lowest 100 mb
  #            raised from 50 mb to 500 mb
  # *** Parcel Theta-E (K) ***
  # *** Lifting Condensation Level (mb) ***
  #

  # Compute parcel surface values
  parcelBoundary  = P[0] - 100.0
  P_parcel  = P[0] - 50.0

  T_sum  = 0.0
  Td_sum  = 0.0
  weightSum = 0.0
  levelTop = 0

  for i in range(1,N):
    if (P[i] > parcelBoundary):
      P_weight  = (P[i-1] - P[i])
      T_sum     = T_sum + (T[i-1] + T[i]) * P_weight/2.0
      Td_sum    = Td_sum + (Td[i-1] + Td[i]) * P_weight/2.0
      weightSum = weightSum + P_weight
      levelTop  = i + 1

  # Compute temperature and dewpoint temperature at parcel top
  T_top = (T[levelTop] + (T[levelTop-1] - T[levelTop]) *
     (parcelBoundary - P[levelTop])/(P[levelTop-1] - P[levelTop]))
  Td_top = (Td[levelTop] + (Td[levelTop-1] - Td[levelTop]) *
     (parcelBoundary - P[levelTop])/(P[levelTop-1] - P[levelTop]))

  # Average these values in with other parcel attributes
  P_weight  = (P[levelTop-1] - parcelBoundary)
  T_sum     = T_sum + (T[levelTop-1] + T_top) * P_weight/2.0
  Td_sum    = Td_sum + (Td[levelTop-1] + Td_top) * P_weight/2.0
  weightSum = weightSum + P_weight

  # Now define parcel attributes including lifted condensation level
  T_parcel    = T_sum/weightSum
  Td_parcel   = Td_sum/weightSum
  lcl         = pressure_at_lcl(P_parcel, T_parcel, Td_parcel, P[0])

  StabilityIndices['LCL']  = lcl
  find_lcl_index   = abs(lcl-P)
  lcl_index  = num.where(find_lcl_index == min(find_lcl_index))

  # Calculate parcel theta-E
  thetaE_parcel = parcel_thetae(P_parcel,T_parcel,Td_parcel,P[0])
  StabilityIndices['Parcel_ThetaE'] = thetaE_parcel[0][0]

  # Calculate lifted index
  StabilityIndices['LiftedIndex'] = (T[N_500] -
             temp_along_satadiabat(thetaE_parcel,500.0,P[0]))[0]

  #
  # *** Showalter Index (degrees C) ***
  # Showalter = T_500 - Tparcel
  # T_parcel  = parcel raised from 850 mb to 500 mb
  #
  thetaE_850 = parcel_thetae(850.0,T[N_850],Td[N_850],P[0])
  StabilityIndices['Showalter'] = (T[N_500] -
      temp_along_satadiabat(thetaE_850,500.0,P[0]))[0]

  #
  # *** Total Totals (degrees C) ***
  # TotalTotals = T_850 + Td_850 - 2*T_500
  #
  StabilityIndices['TotalTotals'] = (T[N_850] + Td[N_850] - 2.0*T[N_500])[0]

  #
  # *** Cross Totals (degrees C) ***
  # CrossTotals = Td_850 - T_500
  #
  StabilityIndices['CrossTotals'] = (Td[N_850] - T[N_500])[0]

  #
  # *** Vertical Totals (degrees C) ***
  # VertTotals = T_850 - T_500
  #
  StabilityIndices['VertTotals'] = (T[N_850] - T[N_500])[0]

  #
  # *** K Index (degrees C) ***
  # KIndex = (T_850 + Td_850) - T_500 - (T_700 - Td_700)
  #
  StabilityIndices['KIndex'] = ((T[N_850] + Td[N_850]) - T[N_500] -
               (T[N_700] - Td[N_700]) - 273.15)[0]

  #
  # *** Level of Free Convection (mb) ***     
  #
  thetae_sounding  = num.zeros((N,))
  tempdiff         = num.zeros((N,))
  N_LFC            = 0
  N_EQ             = 0
  maxthetae_level  = 0

  for i in range(0,N-15):
    thetae_sounding[i] = parcel_thetae(P[i],T[i],Td[i],P[0])

  q = num.where(P > 500)
  maxthetae_level = num.where(thetae_sounding[q] == max(thetae_sounding[q]))
  maxthetae_temp  = thetae_sounding[maxthetae_level]

  for i in range(0,N-15):
    tempdiff[i] = temp_along_satadiabat(thetaE_parcel,P[i],P[0]) - T[i]; 

  for i in range(0,N-15):
    if ((tempdiff[i] > 1.0) and (N_LFC == 0) and (i > maxthetae_level[0])):
      N_LFC = i
    if ((tempdiff[i] < -1.5) and (N_LFC != 0) and (N_EQ == 0)):
      N_EQ = i

  if ((N_LFC != 0) and (N_EQ == 0)):
    N_EQ = N-16

  if ((N_LFC > 0) and (N_LFC < 32)):
    StabilityIndices['LFC'] = P[N_LFC]
  else:
    StabilityIndices['LFC'] = -999.0

  #
  # *** Equilibrium Level (mb) ***
  #
  if ( N_EQ != 0 ):
    StabilityIndices['EQ'] = P[N_EQ]
  else:
    StabilityIndices['EQ'] = -999.0

  #
  # *** CAPE (Convective Available Potential Energy) (J kg^-1) ***
  # *** CIN (Convective Inhibition) (J kg^-1) ***
  #
  [ca,ci] = cape_sound(P,T,R)
  StabilityIndices['CIN']  = ci[0]
  StabilityIndices['CAPE'] = ca[0]

  #
  # *** Total Precipitable Water Vapor (cm) ***
  #
  StabilityIndices['TPWV'] = tpwv(P,R/1000.0)
  return StabilityIndices

if __name__ == '__main__':
  from tdew import tdew
  import pprint
  print(stability.__doc__)
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
  td = tdew(p,t,r)
  Stability = stability(p,t,td)
  pp = pprint.PrettyPrinter(indent=4)
  pp.pprint(Stability)
