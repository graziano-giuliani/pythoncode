#!/usr/bin/env python
# hypsometric.py
import numpy as num

def hypsometric(pressure,temperature,mixing_ratio,surface_altitude):
  """altitude = hypsometric(pressure,temperature,mixing_ratio,surface_altitude)
 
NAME:
      geopotential_altitude   - Paul van Delst's name of the routine
                      (DT renamed it to this).

PURPOSE:
 This function calculates geopotential altitudes using the hypsometric equation.
 It is like the "hypsometric" function I have implemented, but this one
 also accounts for water vapor in the atmosphere.


CATEGORY:
      Meteorology

CALLING SEQUENCE:
      altitude = hypsometric( pressure, $          ; Input
                              temperature, $       ; Input
                              mixing_ratio, $      ; Input
                              surface_altitude, $  ; Input

INPUTS:
      pressure:        Array of pressure in mb
      temperature:     Array of temperature in K
      mixing_ratio:    Array of water vapor mixing ratio in g/kg
      surface_height:  Height of surface in km. 

KEYWORD PARAMETERS:
      None.

OUTPUTS:
      altitude:        Geopotential altitudes in km in the same order as the
                       input data.

CALLS:
      None.


COMMON BLOCKS:
      None

SIDE EFFECTS:
      None

RESTRICTIONS:
   - Input pressure, temperature and mixing_ratio MUST be arrays with at
     LEAST two elements where the element with the highest pressure corresponds
     to the surface altitude passed.
   - Input surface height must be > or = to 0.0km and a SCALAR.
   - No allowance is made, yet, for the change in the acceleration due to 
     gravity with altitude.

PROCEDURE:
      Geopotential heights are calculated using the hypsometric equation:

                        -
                   Rd * Tv    [  p1  ]
        z2 - z1 = --------- ln[ ---- ]
                      g       [  p2  ]

      where Rd    = gas constant for dry air (286.9968933 J/K/kg),
            g     = acceleration due to gravity (9.80616m/s^2),
            Tv    = mean virtual temperature for an atmospheric layer,
            p1,p2 = layer boundary pressures, and
            z1,z2 = layer boundary heights.

      The virtual temperature, the temperature that dry air must have in
      order to have the same density as moist air at the same pressure, is
      calculated using:

                 [      1 - eps      ]
        Tv = T * [ 1 + --------- * w ]
                 [        eps        ]

  where T   = temperature,
        w   = water vapor mixing ratio, and
        eps = ratio of the molecular weights of water and dry air (0.621970585).

EXAMPLE:
      Given arrays of pressure, temperature, and mixing ratio:

        IDL> PRINT, p, t, mr
              1015.42      958.240
              297.180      291.060
              7.83735      5.71762

      the geopotential altitudes can be found by typing:

        IDL> result = geopotential_altitude( p, t, mr, 0.0, alt )
        IDL> PRINT, result, alt
               1      0.00000     0.500970

MODIFICATION HISTORY:
      Written by:     Paul van Delst, CIMSS/SSEC, 08-Dec-1997

Log: hypsometric.pro,m  
converted to matlab code - leslie moy 02-02-02

$Log: hypsometric2.pro,v $
Revision 1.1  2002/01/16 13:28:15  dturner
Initial revision

Revision 1.1  1999/03/03 16:24:55  paulv
Adapted from pressure_height.pro. Improved input argument checking.
  """
  #----------------------------------------------------------------------------
  #                      -- Declare some constants --
  #
  #  Parameter         Description                 Units
  #  ---------   ---------------------------     -----------
  #     Rd       Gas constant for dry air        J/degree/kg
  #
  #     g        Acceleration due to gravity       m/s^2
  #
  #    eps       Ratio of the molec. weights        None
  #              of water and dry air
  #              
  #----------------------------------------------------------------------------
  Rd  = 286.9968933
  g   = 9.80616
  eps = 0.621970585
  #----------------------------------------------------------------------------
  # Calculate average dP
  # --------------------
  n_levels = pressure.size
  dp_average = (num.sum(pressure[0:n_levels-1]-pressure[1:n_levels]) /
               (n_levels))
  # ---------------------------------------------------
  # Sort arrays based on  average pressure differential
  # ---------------------------------------------------
  if ( dp_average < 0.0 ):
    p  = pressure[::-1]
    t  = temperature[::-1]
    mr = mixing_ratio[::-1]
    sort  = 1;
  else:
    p  = pressure
    t  = temperature
    mr = mixing_ratio
    sort  = 0
  #----------------------------------------------------------------------------
  #                   -- Calculate average temperature --
  #----------------------------------------------------------------------------
  t_average = num.zeros((n_levels))
  for i in range(0,n_levels):
    t_average[i] = 0.5 * (t[i] + t[i-1])

  #----------------------------------------------------------------------------
  #       -- Calculate average mixing ratio (in kg/kg - hence --
  #       -- the divisor of 2000.0 instead of 2.0 )           --
  #----------------------------------------------------------------------------
  mr_average = num.zeros((n_levels))
  for i in range(0,n_levels):
    mr_average[i] = 0.0005 * (mr[i] + mr[i-1])

  #----------------------------------------------------------------------------
  #                 -- Calculate virtual temperature --
  #----------------------------------------------------------------------------
  ratio     = (1.0-eps) / eps
  t_virtual = t_average * (1.0+(ratio*mr_average))

  #----------------------------------------------------------------------------
  #         -- Calculate altitudes (divide by 1000.0 to get km) --
  #         -- Make sure the data is returned in an order       --
  #         -- consistent with the input                        --
  #----------------------------------------------------------------------------

  # -- Calculate layer thicknesses
  altitude = num.zeros((n_levels))
  altitude[0] = surface_altitude
  for i in range(1,n_levels):
    altitude[i] = (altitude[i-1] + 0.001 * 
           ((Rd/g)*t_virtual[i]*num.log(p[i-1]/p[i])))

  # -- Reverse array order if required
  if ( sort == 1 ):
    altitude = altitude[::-1]

  return altitude

if __name__ == '__main__':
  print(hypsometric.__doc__)
  p = num.array(
      ( 1012.0, 991.3, 969.1, 945.5, 920.4, 893.8, 865.7, 836.1, 805.1, 772.8,
         739.5, 705.2, 670.3, 635.0, 599.7, 564.5, 529.8, 495.7, 462.6, 430.7,
         400.0, 370.8, 343.0, 316.7, 292.0, 266.8, 247.2, 227.0, 208.2, 190.8,
         174.7, 159.9, 146.2, 133.6, 121.9, 111.3, 101.5,  92.6,  84.4,  76.9,
          70.0 ))
  t = num.array(
      ( 24.54, 23.16, 21.67, 20.23, 18.86, 17.49, 16.10, 14.69, 13.22, 11.52,
         9.53,  7.24,  4.80,  2.34,  0.04, -2.29, -4.84, -7.64,-10.66,-13.95,
       -17.54,-21.45,-25.58,-29.90,-34.33,-38.94,-43.78,-48.80,-53.94,-58.79,
       -63.27,-67.32,-70.74,-73.62,-75.74,-77.07,-77.43,-76.63,-75.06,-73.14,
       -71.43 ))
  t = t + 273.15
  r = num.array(
      ( 17.78, 16.92, 15.93, 14.87, 13.78, 12.70, 11.84, 10.96, 10.15,  9.31,
         8.46,  7.73,  7.05,  6.32,  5.62,  4.91,  4.10,  3.30,  2.67,  2.15,
         1.66,  1.26,  0.95,  0.68,  0.45,  0.28,  0.17,  0.10,  0.06,  0.04,
         0.02,  0.02,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,
         0.02 ))
  print(p)
  print(t)
  print(r)
  h = hypsometric(p,t,r,0.001)
  print(h)
  p = p[::-1]
  t = t[::-1]
  r = r[::-1]
  h = hypsometric(p,t,r,0.001)
  print(h)
