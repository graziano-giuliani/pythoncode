#!/usr/bin/env python
# qtips.py
import numpy as num

def qtips(gid):
  """[a,b,c,d,g] = qtips(gid)

returns nuclear degeneracy factors g and Gamache's internal partition
sum coefficients a, c, c, and d.  This is taken from Dave Edwards 
qtips.f.

***********************************************************************

  PROGRAM        QTIPS   BLOCK DATA

  PURPOSE        TOTAL INTERNAL PARTITION SUMS

  VERSION        3.0   D.P. EDWARDS   11/01/90 

  DESCRIPTION    COEFFICIENTS FOR CALCULATING THE TOTAL INTERNAL 
                 PARTITION SUMS
  
                 THIS IS A COPY OF THE ROUTINE QSUMS
                 FROM THE PROGRAM TIPS BY R.R.GAMACHE WHICH CALCULATES
                 THE TOTAL INTERNAL PARTITION FUNCTIONS IN THE 
                 HITRAN SELECT ROUTINE
                 ...LAST MODIFIED MAY 24, 1991

***********************************************************************

...STATE INDEPENDENT DEGENERACY FACTORS: GJ 
...(INCLUDES GENERAL NUCLEAR FACTOR (P(2I+1)), (2S+1), AND (2-DL0)

     DATA GJ/
  """
  amus_mols = ( 'H2O',  'CO2',  'O3',   'N2O',  'CO',   'CH4', 'O2',  'NO',
                'SO2',  'NO2',  'NH3',  'HNO3', 'OH',   'HF',  'HCL', 'HBR',
                'HI',   'CLO',  'OCS',  'H2CO', 'HOCL', 'N2',  'HCN', 'CH3CL',
                'H2O2', 'C2H2', 'C2H6', 'PH3',  'COF2', 'SF6', 'H2S', 'HCOOH' )
  try:
    gasid = amus_mols.index(gid)+1
  except:
    gasid = gid
  if (gasid ==  1):
    g = num.array(( 1, 1, 6, 6 )) # H2O
  elif (gasid ==  2):
    g = num.array(( 1, 2, 1, 6, 2, 12, 1, 6 ))  # CO2
  elif (gasid ==  3):
    g = num.array(( 1, 1, 1 )) # O3
  elif (gasid ==  4):
    g = num.array(( 9, 6, 6, 9,54 )) # N2O
  elif (gasid ==  5):
    g = num.array(( 1, 2, 1, 6,2 )) # CO
  elif (gasid ==  6):
    g = num.array(( 1, 2, 3 )) # CH4
  elif (gasid ==  7):
    g = num.array(( 1, 1, 6 )) # O2  
  elif (gasid ==  8):
    g = num.array(( 12, 8, 12 )) # NO
  elif (gasid ==  9):
    g = num.array(( 1, 1 )) # SO2
  elif (gasid == 10):
    g = num.array(( 6 )) # NO2
  elif (gasid == 11):
    g = num.array(( 3, 2 )) # NH3
  elif (gasid == 12):
    g = num.array(( 6 )) # HNO3
  elif (gasid == 13):
    g = num.array(( 8, 8, 12 )) # OH
  elif (gasid == 14):
    g = num.array(( 4 )) # HF
  elif (gasid == 15):
    g = num.array(( 8, 8 )) # HCL
  elif (gasid == 16):
    g = num.array(( 8, 8 )) # HBR
  elif (gasid == 17):
    g = num.array(( 12 )) # HI
  elif (gasid == 18):
    g = num.array(( 4, 4 )) # CLO
  elif (gasid == 19):
    g = num.array(( 1, 1, 2, 1 )) # OCS
  elif (gasid == 20):
    g = num.array(( 1, 2, 1 )) # H2CO
  elif (gasid == 21):
    g = num.array(( 8, 8 )) # HOCL
  elif (gasid == 22):
    g = num.array(( 0.5 )) # N2 
  elif (gasid == 23):
    g = num.array(( 6, 12, 4 )) # HCN
  elif (gasid == 24):
    g = num.array(( 4, 4 )) # CH3CL
  elif (gasid == 25):
    g = num.array(( 4 )) # H2O2
  elif (gasid == 26):
    g = num.array(( 1, 8 )) # C2H2
  elif (gasid == 27):
    g = num.array(( 64 )) # C2H6
  elif (gasid == 28):
    g = num.array(( 2 )) # PH3
  elif (gasid == 29):
    g = num.array(( 1 )) # COF2
  elif (gasid == 30):
    g = num.array(( 1 )) # SF6
  elif (gasid == 31):
    g = num.array(( 1 )) # H2S
  elif (gasid == 32):
    g = num.array(( 1 )) # HCOOH
  else:
    raise ValueError('gasid not in range or unkown')

  g = num.reshape(g,(-1,1))

  #...TOTAL INTERNAL PARTITION SUMS FOR 70 - 405 K RANGE

  if (gasid == 1):    # isotopes 161  181  171  162
    abcd = num.array((
      -.37688E+01, .26168E+00, .13497E-02, -.66013E-06,
      -.38381E+01, .26466E+00, .13555E-02, -.65372E-06,
      -.22842E+02, .15840E+01, .81575E-02, -.39650E-05,
      -.20481E+02, .13017E+01, .66225E-02, -.30447E-05))
  elif (gasid == 2):
    abcd = num.array((
      -.21995E+01, .96751E+00, -.80827E-03, .28040E-05,
      -.38840E+01, .19263E+01, -.16058E-02, .58202E-05,
      -.47289E+01, .20527E+01, -.17421E-02, .60748E-05,
      -.27475E+02, .11973E+02, -.10110E-01, .35187E-04,
      -.84191E+01, .41186E+01, -.34961E-02, .12750E-04,
      -.48468E+02, .23838E+02, -.20089E-01, .73067E-04,
      -.22278E+01, .10840E+01, -.89718E-03, .32143E-05,
      -.29547E+02, .12714E+02, -.10913E-01, .38169E-04))
  elif (gasid == 3):
    abcd = num.array((
      -.13459E+03, .62255E+01,  .14811E-01, .18608E-04,
      -.12361E+03, .61656E+01,  .19168E-01, .13223E-04,
      -.12359E+03, .60957E+01,  .18239E-01, .13939E-04))
  elif (gasid == 4):
    abcd = num.array((
      -.95291E+01, .15719E+02, -.12063E-01, .53781E-04,
       .48994E+01, .10211E+02, -.62964E-02, .33355E-04,
      -.28797E+01, .10763E+02, -.78058E-02, .36321E-04,
       .25668E+02, .15803E+02, -.67882E-02, .44093E-04,
       .18836E+03, .91152E+02, -.31071E-01, .23789E-03))
  elif (gasid == 5):
    abcd = num.array((
       .31591E+00, .36205E+00, -.22603E-05, .61215E-08,
       .62120E+00, .75758E+00, -.59190E-05, .15232E-07,
       .30985E+00, .38025E+00, -.29998E-05, .76646E-08,
       .18757E+01, .22289E+01, -.15793E-04, .41607E-07,
       .60693E+00, .79754E+00, -.78021E-05, .19200E-07))
  elif (gasid == 6):
    abcd = num.array((
      -.17475E+02, .95375E+00,  .39758E-02,-.81837E-06,
      -.27757E+02, .17264E+01,  .93304E-02,-.48181E-05,
      -.89810E+03, .44451E+02,  .17474E+00,-.22469E-04))
  elif (gasid == 7):
    abcd = num.array((
      -.10000E+01, .00000E+00,  .00000E+00, .00000E+00,
      -.10000E+01, .00000E+00,  .00000E+00, .00000E+00,
      -.10000E+01, .00000E+00,  .00000E+00, .00000E+00))
  elif (gasid == 8):
    abcd = num.array((
      -.17685E+03, .28839E+02,  .87413E-01,-.92142E-04,
      -.61157E+02, .13304E+02,  .40161E-01,-.42247E-04,
      -.18775E+03, .30428E+02,  .92040E-01,-.96827E-04))
  elif (gasid == 9):
    abcd = num.array((
      -.17187E+03, .94104E+01,  .34620E-01, .25199E-04,
      -.17263E+03, .94528E+01,  .34777E-01, .25262E-04))
  elif (gasid == 10):
    abcd = num.array((
      -.89749E+03, .44718E+02,  .15781E+00, .43820E-04))
  elif (gasid == 11):
    abcd = num.array((
      -.48197E+02, .27739E+01,  .11492E-01,-.18209E-05,
      -.32700E+02, .18444E+01,  .77001E-02,-.12388E-05))
  elif (gasid == 12):
    abcd = num.array((
      -.74208E+04, .34984E+03,  .89051E-01, .39356E-02))
  elif (gasid == 13):
    abcd = num.array((
       .76510E+02, .11377E+01,  .39068E-02,-.42750E-05,
       .76140E+02, .11508E+01,  .39178E-02,-.42870E-05,
       .14493E+03, .47809E+01,  .15441E-01,-.16217E-04))
  elif (gasid == 14):
    abcd = num.array((
       .15649E+01, .13318E+00,  .80622E-05,-.83354E-08))
  elif (gasid == 15):
    abcd = num.array((
       .28877E+01, .53077E+00,  .99904E-05,-.70856E-08,
       .28873E+01, .53157E+00,  .99796E-05,-.70647E-08))
  elif (gasid == 16):
    abcd = num.array((
       .28329E+01, .66462E+00,  .83420E-05,-.30996E-08,
       .28329E+01, .66483E+00,  .83457E-05,-.31074E-08))
  elif (gasid == 17):
    abcd = num.array((
       .41379E+01, .12977E+01,  .61598E-05, .10382E-07))
  elif (gasid == 18):
    abcd = num.array((
       .15496E+04, .11200E+03,  .19225E+00, .40831E-04,
       .15728E+04, .11393E+03,  .19518E+00, .43308E-04))
  elif (gasid == 19):
    abcd = num.array((
       .18600E+02, .31185E+01,  .30405E-03, .85400E-05,
       .19065E+02, .31965E+01,  .31228E-03, .87535E-05,
       .42369E+02, .61394E+01,  .13090E-02, .16856E-04,
       .21643E+02, .32816E+01,  .57748E-03, .90034E-05))
  elif (gasid == 20):
    abcd = num.array((
      -.44663E+02, .23031E+01,  .95095E-02,-.16965E-05,
      -.91605E+02, .47223E+01,  .19505E-01,-.34832E-05,
      -.44663E+02, .23031E+01,  .95095E-02,-.16965E-05))
  elif (gasid == 21):
    abcd = num.array((
      -.62547E+03, .31546E+02,  .11132E+00, .32438E-04,
      -.60170E+03, .31312E+02,  .11841E+00, .23717E-04))
  elif (gasid == 22):
    abcd = num.array((
       .73548E+00, .78662E+00, -.18282E-05, .68772E-08))
  elif (gasid == 23):
    abcd = num.array((
      -.97107E+00, .29506E+01, -.16077E-02, .61148E-05,
      -.16460E+01, .60490E+01, -.32724E-02, .12632E-04,
      -.40184E+00, .20202E+01, -.10855E-02, .42504E-05))
  elif (gasid == 24):
    abcd = num.array((
      -.89695E+03, .40155E+02,  .82775E-01, .13400E-03,
      -.91113E+03, .40791E+02,  .84091E-01, .13611E-03))
  elif (gasid == 25):
    abcd = num.array((
      -.95255E+03, .49483E+02,  .21249E+00,-.35489E-04))
  elif (gasid == 26):
    abcd = num.array((
       .25863E+01, .11921E+01, -.79281E-03, .46225E-05,
       .20722E+02, .95361E+01, -.63398E-02, .36976E-04))
  elif (gasid == 27):
    abcd = num.array((
      -.10000E+01, .00000E+00,  .00000E+00, .00000E+00))
  elif (gasid == 28):
    abcd = num.array((
      -.11388E+03, .69602E+01,  .17396E-01, .65088E-05))
  elif (gasid == 29):
    abcd = num.array((
      -.10000E+01, .00000E+00,  .00000E+00, .00000E+00))
  elif (gasid == 30):
    abcd = num.array((
      -.10000E+01, .00000E+00,  .00000E+00, .00000E+00))
  elif (gasid == 31):
    abcd = num.array((
      -.10000E+01, .00000E+00,  .00000E+00, .00000E+00))
  elif (gasid == 32):
    abcd = num.array((
      -.10000E+01, .00000E+00,  .00000E+00, .00000E+00))
  else:
    raise ValueError('gasid not in range or unkown')
  a = num.reshape(abcd[0::4],(-1,1))
  b = num.reshape(abcd[1::4],(-1,1))
  c = num.reshape(abcd[2::4],(-1,1))
  d = num.reshape(abcd[3::4],(-1,1))
  return [a,b,c,d,g]

if __name__ == '__main__':
  print(qtips.__doc__)
  [a,b,c,d,g] = qtips(1)
  print(a)
  print(b)
  print(c)
  print(d)
  print(g)
  try:
    mm = gas_mass(-14)
  except:
    print('Ok')
  [a,b,c,d,g] = qtips('CO2')
  print(a)
  print(b)
  print(c)
  print(d)
  print(g)
