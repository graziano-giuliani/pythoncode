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
import numpy as np

class oss4SHIS:
  """Use HITRAN precomputed tables to compute oss forward model (ir)"""
  def __init__(self,solar,od,imolid=None,imolind=None):
    from asolar import asolar
    from hitran import hitran
    from ossir import oss_ir
    a = asolar(solar)
    o = hitran(od)
    self.oss = oss_ir
    if ( imolid is not None ):
      uimolid = imolid
    else:
      uimolid = np.array([1,2,3],np.int32)
    if ( imolind is not None ):
      self.oss.set_imols(uimolid,imolind)
    else:
      self.oss.set_imols(uimolid)
    pref = o.get('pref')
    tmptab = o.get('tmptab')
    self.oss.set_hitran(pref,tmptab)
    self.cwvn = o.get('cWvn')
    vwvn = o.get('vwvn')
    sunrad = a.get(vwvn)
    self.oss.set_solar_irradiance(vwvn,sunrad)
    wvptab = o.get('wvptab')
    kfix = o.get('kfix')
    dkh2o = o.get('dkh2o')
    kh2o = o.get('kh2o')
    kvar = o.get('kvar')
    imols = o.get('imols')
    molid = o.get('molid')
    nch = o.get('nch')
    coef = o.get('coef')
    ichmap = o.get('ichmap')
    self.nfs = o.nf()
    self.nchan = o.nchan()
    self.oss.set_hitran_absorption_coefficients(molid,imols,
        wvptab,coef,kfix,kh2o,dkh2o,kvar,nch,ichmap)
  def compute(self,indata,outdata):
    if ( 'xg' in indata.keys() ):
      xg = indata['xg']
      initempg = indata['initempg']
      initsking = indata['initsking']
      inipsfcg = indata['inipsfcg']
    elif ( 'temp' in indata.keys() ):
      nlev = len(indata['temp'])
      xg = np.zeros(nlev*4+2, dtype=float)
      xg[0:nlev] = indata['temp']
      xg[nlev] = indata['tskin']
      xg[nlev+1] = indata['psf']
      xg[nlev+2:nlev+2+nlev] = indata['h2o']
      if ( 'co2' not in indata.keys() ):
        xg[2*nlev+2:2*nlev+2+nlev] = ( 
           np.zeros(nlev,dtype=np.float32,order='F') + 5.9407297521829605E-004 )
      else:
        xg[2*nlev+2:2*nlev+2+nlev] = indata['co2']
      xg[3*nlev+2:3*nlev+2+nlev] = indata['o3']
      # Fortran indexes below !
      initempg  = 1
      initsking = nlev+1
      inipsfcg  = initsking+1
    else:
      raise KeyError('Neither xg or temperature data in input')
    pobs = indata['pobs']
    obsang = indata['obsang']
    sunang = indata['sunang']
    if ( 'y' not in outdata.keys() ):
      nparg = len(xg)
      outdata['y'] = np.zeros(self.nchan,np.float32,order='F')
      outdata['xkt'] = np.zeros((nparg,self.nchan),np.float32,order='F')
      outdata['xkemrf'] = np.zeros((2,self.nfs,self.nchan),np.float32,order='F')
      outdata['paxkemrf'] = np.zeros((2,self.nchan),np.float32,order='F')
    if ( 'pressure' in indata.keys() ):
      self.oss.ossdrv_ir(initempg,initsking,inipsfcg,xg,pobs,obsang,sunang,
         indata['sfgrd'],indata['emrf'],
         outdata['y'],outdata['xkt'],outdata['xkemrf'],outdata['paxkemrf'],
         puser=indata['pressure'])
    else:
      self.oss.ossdrv_ir(initempg,initsking,inipsfcg,xg,pobs,obsang,sunang,
         indata['sfgrd'],indata['emrf'],
         outdata['y'],outdata['xkt'],outdata['xkemrf'],outdata['paxkemrf'])
  def __del__(self):
    self.oss.release( )

#
# Unit test of the above class
#
if ( __name__ == '__main__' ):
  from netCDF4 import Dataset
  from emissivity import emissivity

  # Initialize oss forward model
  a = oss4SHIS('../data/solar_irradiances.nc',
               '../data/leo.iasi.0.05.nc')

  # Input data
  e = emissivity('../data/emissivity.nc')
  indata = {}
  indata['sfgrd'] = e.get('SfGrd')
  indata['emrf'] = e.get('EmRf')
  indata['tskin'] = 299.82998657226562
  indata['psf'] = 997.13000488281250
  indata['temp'] = np.array([
   195.49699401855469 , 205.28900146484375 , 215.74800109863281 ,
   229.22200012207031 , 241.72200012207031 , 252.51400756835938 ,
   262.03500366210938 , 269.09399414062500 , 270.45498657226562 ,
   268.20300292968750 , 263.84600830078125 , 258.26800537109375 ,
   255.75300598144531 , 251.98300170898438 , 243.87500000000000 ,
   238.94999694824219 , 236.05299377441406 , 233.47799682617188 ,
   229.43099975585938 , 226.70899963378906 , 225.07699584960938 ,
   222.82899475097656 , 220.17999267578125 , 218.31599426269531 ,
   217.51499938964844 , 216.85499572753906 , 216.33599853515625 ,
   216.05499267578125 , 216.19000244140625 , 215.85600280761719 ,
   215.19999694824219 , 214.99699401855469 , 214.36900329589844 ,
   213.96699523925781 , 214.48399353027344 , 215.22399902343750 ,
   215.52400207519531 , 215.58399963378906 , 215.48599243164062 ,
   214.65600585937500 , 214.46400451660156 , 214.98800659179688 ,
   214.32099914550781 , 213.68600463867188 , 214.66700744628906 ,
   216.43600463867188 , 216.22099304199219 , 214.95899963378906 ,
   214.44599914550781 , 214.70899963378906 , 215.75100708007812 ,
   217.52999877929688 , 219.94099426269531 , 222.82200622558594 ,
   225.99099731445312 , 229.28300476074219 , 232.70799255371094 ,
   236.20599365234375 , 239.64300537109375 , 243.25300598144531 ,
   247.07200622558594 , 250.92399597167969 , 254.60299682617188 ,
   257.96398925781250 , 260.90100097656250 , 263.36898803710938 ,
   265.57598876953125 , 267.63101196289062 , 269.60101318359375 ,
   271.66101074218750 , 273.79098510742188 , 275.86599731445312 ,
   277.73901367187500 , 279.39199829101562 , 280.90399169921875 ,
   282.29299926757812 , 283.55999755859375 , 284.77301025390625 ,
   285.93499755859375 , 286.76199340820312 , 287.62298583984375 ,
   288.86199951171875 , 290.03799438476562 , 291.07699584960938 ,
   291.98001098632812 , 292.75399780273438 , 293.42300415039062 ,
   294.00601196289062 , 294.49899291992188 , 294.92199707031250 ,
   295.37899780273438 ] , dtype=np.float32,order='F')
  indata['h2o'] = np.array([
   1.8522099480833276E-006 , 2.5321699013147736E-006 , 3.3462299597886158E-006 ,
   4.0941399674920831E-006 , 4.1902098928403575E-006 , 4.1866401261358988E-006 ,
   4.1710000004968606E-006 , 4.1224798223993275E-006 , 4.0048498703981750E-006 ,
   3.8550301724171732E-006 , 3.7657900975318626E-006 , 3.6219200865161838E-006 ,
   3.4851100281230174E-006 , 3.4485299238440348E-006 , 3.4380600482109003E-006 ,
   3.4100401080650045E-006 , 3.3074200018745614E-006 , 3.1797599149285816E-006 ,
   3.0800499644101365E-006 , 3.0216399409255246E-006 , 2.9836101020919159E-006 ,
   2.9567499950644560E-006 , 2.9416301003948320E-006 , 2.9327600259421160E-006 ,
   2.9243399239931023E-006 , 2.9125098990334664E-006 , 2.8948700219189050E-006 ,
   2.8727899916702881E-006 , 2.8523600121843629E-006 , 2.8193001071485924E-006 ,
   2.7724099709303118E-006 , 2.7197600047657033E-006 , 2.6764701033243909E-006 ,
   2.6500499643589137E-006 , 2.6541599709162256E-006 , 2.6728100692707812E-006 ,
   2.6880099994741613E-006 , 2.7176899948244682E-006 , 2.7499499992700294E-006 ,
   2.7537000732991146E-006 , 2.9873099265387282E-006 , 3.4484698971937178E-006 ,
   4.3864301915164106E-006 , 5.7857801039062906E-006 , 7.8530601967941038E-006 ,
   9.5091900220722891E-006 , 1.0954900062642992E-005 , 1.3077799849270377E-005 ,
   1.6240199329331517E-005 , 2.1141200704732910E-005 , 2.9084299967507832E-005 ,
   4.0059898310573772E-005 , 5.5496500863227993E-005 , 7.7732700447086245E-005 ,
   1.0551000013947487E-004 , 1.3935800234321505E-004 , 1.8844100122805685E-004 ,
   2.4947800557129085E-004 , 2.8610401204787195E-004 , 2.5873200502246618E-004 ,
   2.2800000442657620E-004 , 2.2262299899011850E-004 , 2.3247400531545281E-004 ,
   2.6885300758294761E-004 , 3.6187001387588680E-004 , 5.0088198622688651E-004 ,
   6.4087699865922332E-004 , 8.5206201765686274E-004 , 1.2005500029772520E-003 ,
   1.7163499724119902E-003 , 2.2430599201470613E-003 , 2.7267099358141422E-003 ,
   3.3067099284380674E-003 , 4.0873899124562740E-003 , 4.9571897834539413E-003 ,
   5.8562001213431358E-003 , 6.5321098081767559E-003 , 6.8273600190877914E-003 ,
   6.9653899408876896E-003 , 7.3573398403823376E-003 , 8.1664798781275749E-003 ,
   8.4378998726606369E-003 , 8.4838503971695900E-003 , 8.5134999826550484E-003 ,
   8.5416501387953758E-003 , 8.5704904049634933E-003 , 8.6027402430772781E-003 ,
   8.6434800177812576E-003 , 8.7000103667378426E-003 , 8.7922401726245880E-003 ,
   9.0125501155853271E-003 ] , dtype=np.float32,order='F')
  indata['o3'] = np.array([
   2.0679999579442665E-007 , 3.7505398609027907E-007 , 7.1199701778823510E-007 ,
   1.2760500567310373E-006 , 1.8921000446425751E-006 , 2.4401199425483355E-006 ,
   2.9400000585155794E-006 , 3.4252600471518235E-006 , 4.0483300836058334E-006 ,
   4.9532200137036853E-006 , 6.1584300965478178E-006 , 7.7594104368472472E-006 ,
   9.3812495833844878E-006 , 1.1037899639632087E-005 , 1.2533399967651349E-005 ,
   1.3236700397101231E-005 , 1.3621300240629353E-005 , 1.3695699635718483E-005 ,
   1.3193500308261719E-005 , 1.2132600204495247E-005 , 1.0765000297396909E-005 ,
   9.8807104222942144E-006 , 9.1616902864188887E-006 , 8.5475803643930703E-006 ,
   7.8542198025388643E-006 , 7.2931002250697929E-006 , 6.7923201640951447E-006 ,
   6.1678101701545529E-006 , 5.5271798373723868E-006 , 4.9240297812502831E-006 ,
   4.2081301216967404E-006 , 3.2737400488258572E-006 , 2.4478299565089401E-006 ,
   1.8292099639438675E-006 , 1.5723200021966477E-006 , 1.7079499912142637E-006 ,
   1.8751200059341500E-006 , 1.8040799432128551E-006 , 1.4608200444854447E-006 ,
   1.2204500308143906E-006 , 1.0136600394616835E-006 , 7.5832298307432211E-007 ,
   4.8743902425485430E-007 , 3.0938198847252352E-007 , 2.8449500177885056E-007 ,
   2.8604600288417714E-007 , 3.2440200925520912E-007 , 3.1362500862996967E-007 ,
   2.2775900276883476E-007 , 1.5017999999145104E-007 , 1.0815799811325633E-007 ,
   9.6226798973475525E-008 , 8.4470997308017104E-008 , 7.6241796875820000E-008 ,
   7.2805399042863428E-008 , 7.0336696467165893E-008 , 6.8959401744450588E-008 ,
   7.1598599049593759E-008 , 7.8118098656432267E-008 , 7.7721900026972435E-008 ,
   7.1205796814410860E-008 , 6.7037298379091226E-008 , 6.4800296684097702E-008 ,
   6.4714200220805651E-008 , 6.5708597674074554E-008 , 6.6245696928035613E-008 ,
   6.6256198749670148E-008 , 6.6784402008579491E-008 , 6.7629699174176494E-008 ,
   6.8265897823494015E-008 , 6.8371399208899675E-008 , 6.8178302115029510E-008 ,
   6.7982902862695482E-008 , 6.8169299538567429E-008 , 6.8325398672186566E-008 ,
   6.8409001130476099E-008 , 6.9107201738916046E-008 , 6.9786501910584775E-008 ,
   7.0195596890698653E-008 , 7.0625702619508957E-008 , 7.1145201729905239E-008 ,
   7.1606102380883385E-008 , 7.2118901073281449E-008 , 7.2543699047855625E-008 ,
   7.2642002635348035E-008 , 7.2760300895424734E-008 , 7.3021503510517505E-008 ,
   7.3266797073756607E-008 , 7.3329900374119461E-008 , 7.3690600288500718E-008 ,
   7.4984598086302867E-008 ] , dtype=np.float32,order='F')
  indata['pressure'] = np.array([
  1.0000200010836124E-002 , 2.9904400929808617E-002 , 5.6840099394321442E-002 ,
  0.10147800296545029 , 0.17160999774932861 , 0.27683201432228088 ,
  0.42849698662757874 , 0.63957101106643677 , 0.92441600561141968 ,
   1.2985099554061890 , 1.7781200408935547 , 2.3799700736999512 ,
   3.1208999156951904 , 4.0175499916076660 , 5.0860300064086914 ,
   6.3416600227355957 , 7.7987999916076660 , 9.4705600738525391 ,
   11.368800163269043 , 13.503700256347656 , 15.884400367736816 ,
   18.517900466918945 , 21.410200119018555 , 24.565299987792969 ,
   27.986000061035156 , 31.673599243164062 , 35.628101348876953 ,
   39.848098754882812 , 44.331001281738281 , 49.073200225830078 ,
   54.070098876953125 , 59.314998626708984 , 64.797798156738281 ,
   70.506103515625000 , 76.429100036621094 , 82.571998596191406 ,
   88.957496643066406 , 95.614097595214844 , 102.57499694824219 ,
   109.87999725341797 , 117.57399749755859 , 125.71199798583984 ,
   134.34599304199219 , 143.51199340820312 , 153.23899841308594 ,
   163.55700683593750 , 174.49699401855469 , 186.09199523925781 ,
   198.37600708007812 , 211.38600158691406 , 225.15800476074219 ,
   239.73300170898438 , 255.15199279785156 , 271.45800781250000 ,
   288.69799804687500 , 306.91699218750000 , 326.16598510742188 ,
   346.49600219726562 , 367.96099853515625 , 390.61700439453125 ,
   414.51599121093750 , 439.69299316406250 , 466.04699707031250 ,
   493.34100341796875 , 521.32000732421875 , 549.71197509765625 ,
   578.24102783203125 , 606.68798828125000 , 634.90899658203125 ,
   662.77899169921875 , 690.18402099609375 , 716.98498535156250 ,
   743.04101562500000 , 768.23797607421875 , 792.47399902343750 ,
   815.63500976562500 , 837.59997558593750 , 858.28399658203125 ,
   877.62402343750000 , 895.54699707031250 , 911.97900390625000 ,
   926.88299560546875 , 940.25402832031250 , 952.09002685546875 ,
   962.38500976562500 , 971.16900634765625 , 978.69702148437500 ,
   985.11401367187500 , 990.30999755859375 , 994.27001953125000 ,
   997.13000488281250], dtype=np.float32,order='F')
  indata['pobs'] = 1.0000200010836124E-002
  indata['obsang'] = 43.695999145507812
  indata['sunang'] = 89.0

  # Output data (empty first call)
  outdata = {}

  #for i in range(0,200):
  #  print(i)
  #  a.compute(indata,outdata)
  a.compute(indata,outdata)

  nlev = len(indata['temp'])
  rootgrp = Dataset('test.nc', 'w', format='NETCDF4')
  nchand = rootgrp.createDimension('channels', a.nchan)
  nsped = rootgrp.createDimension('nspe', a.nfs)
  npargd = rootgrp.createDimension('levels', nlev)
  jj = rootgrp.createDimension('jj', 2)
  jj = rootgrp.createDimension('sfc', 1)
  ncx = rootgrp.createVariable('x','f4',('channels',))
  ncy = rootgrp.createVariable('y','f4',('channels',))
  nckt = rootgrp.createVariable('Kt','f4',('channels','levels'))
  nckwv = rootgrp.createVariable('Kwv','f4',('channels','levels'))
  nckco2 = rootgrp.createVariable('Kco2','f4',('channels','levels'))
  ncko3 = rootgrp.createVariable('Ko3','f4',('channels','levels'))
  nckskt = rootgrp.createVariable('Kskt','f4',('channels','sfc'))
  ncksp = rootgrp.createVariable('Ksp','f4',('channels','sfc'))
  ncxkemrf = rootgrp.createVariable('xkemrf','f4',('channels','nspe','jj',))
  ncpaxkemrf = rootgrp.createVariable('paxkemrf','f4',('channels','jj',))
  ncx[:] = a.cwvn
  ncy[:] = outdata['y']
  nckt[:] = np.transpose(outdata['xkt'][0:nlev,:])
  nckskt[:] = np.transpose(outdata['xkt'][nlev,:])
  ncksp[:] = np.transpose(outdata['xkt'][nlev+1,:])
  nckwv[:] = np.transpose(outdata['xkt'][nlev+2:nlev+2+nlev,:])
  nckco2[:] = np.transpose(outdata['xkt'][2*nlev+2:2*nlev+2+nlev,:])
  ncko3[:] = np.transpose(outdata['xkt'][3*nlev+2:3*nlev+2+nlev,:])
  ncxkemrf[:] = np.transpose(outdata['xkemrf'])
  ncpaxkemrf[:] = np.transpose(outdata['paxkemrf'])
  rootgrp.close()
