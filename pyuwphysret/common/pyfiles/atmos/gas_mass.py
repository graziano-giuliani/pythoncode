#!/usr/bin/env python
# gas_mass.py
from define_constants import *

def gas_mass(gasid):
  """mass_molecule = gas_mass(gasid)

gas_mass.m returns the mass of the HITRAN gas ID gasid

DCT 3/2/1996
  """
  # 
  # note: results are not accurate because amu 
  # values need more significant figures
  # 

  amus_mols = ( 'H2O',  'CO2',  'O3',   'N2O',  'CO',   'CH4', 'O2',  'NO',
                'SO2',  'NO2',  'NH3',  'HNO3', 'OH',   'HF',  'HCL', 'HBR',
                'HI',   'CLO',  'OCS',  'H2CO', 'HOCL', 'N2',  'HCN', 'CH3CL',
                'H2O2', 'C2H2', 'C2H6', 'PH3',  'COF2', 'SF6', 'H2S', 'HCOOH',
                'DAIR' )
  amus_list = ( 2*1 + 16      , 12 + 2*16    , 3*16        , 2*14 + 16       ,
                12 + 16       , 12 + 4*1     , 2*16        , 14 + 16         ,
                32 + 2*16     , 14 + 2*16    , 14 + 3*1    , 1 + 14 + 3*16   ,
                16 + 1        , 1 + 19       , 1 + 35      , 1 + 80          ,
                1 + 127       , 35 + 16      , 16 + 12 + 32, 2*1 + 12 + 16   ,
                1 + 16 + 35   , 2*14         , 1 + 12 + 14 , 12 + 3*1 + 35   ,
                2*1 + 2*16    , 2*12 + 2*1   , 2*12 + 6*1  , 31 + 3*1        ,
                12 + 16 + 2*19, 32 + 6*19    , 2*1 + 32    , 1 + 12 + 2*16 +1,
                28.9402753668809 ) # air (makes water/air=0.621970585)
  try:
    amus = amus_list[gasid-1]
  except:
    if ( gasid == 99 ):
      amus = amus_list[-1]
    else:
      try:
        amus = amus_list[amus_mols.index(gasid)]
      except:
        raise ValueError('gasid not in range or unkown')
  mass_molecule = mass_proton*amus
  return mass_molecule

if __name__ == '__main__':
  print(gas_mass.__doc__)
  mm = gas_mass(4)
  print(mm)
  mm = gas_mass(99)
  print(mm)
  try:
    mm = gas_mass(-14)
  except:
    print('Ok')
  mm = gas_mass('CO2')
  print(mm)
