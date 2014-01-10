#!/usr/bin/env python
# mirto_code_main.py
import numpy as np
from scipy.linalg import lu , solve
import cProfile, pstats
import mirto_code_configuration
import mirto_code_compute_F
import sys

# Empty class to store result
class mirto_state:
  pass

class mirto_residuals:
  pass

class mirto_history:
  def __init__(self):
    self.measurement_space_only = []

class mirto_norm:
  def __init__(self):
    self.measurement_space_only = np.NAN

class mirto:
  def __init__(self,datapath,oss):
    # Load configuration parameters and Input data for retrieval
    self.cx = mirto_code_configuration.mirto_config(datapath,oss)
    self.obsErr = mirto_code_configuration.mirto_obsErr(self.cx)
    self.apriori = mirto_code_configuration.mirto_apriori(self.cx)
    self.state = mirto_state()
    self.norm = mirto_norm()
    self.hist = mirto_history()
    self.residuals = mirto_residuals()
    self.converge = 0.03 * len(self.cx.state_var_indx)

  def compute_chi_square(self):
    self.norm.measurement_space_only = (
      np.dot(np.dot(self.residuals.yobs_minus_yhat.T,self.obsErr.SeInv),
          self.residuals.yobs_minus_yhat/len(self.residuals.yobs_minus_yhat)) )

  def update_solution(self,fm):
    KtSeInv = np.dot(fm.K.T,self.obsErr.SeInv)
    KtSeInvK = np.dot(KtSeInv,fm.K)
    A = (KtSeInvK+(1.0+self.cx.gamma)*self.state.SaInv_ret)
    dx = (self.state.xhat-self.state.xa)
    d = (np.dot(KtSeInv,self.residuals.yobs_minus_yhat) - 
         np.dot(self.state.SaInv_ret,dx))
    # Use iterative LU decomposition to determine the solution
    # First iteration
    L,U = lu(A,permute_l=True)
    y = solve(L,d)
    x = solve(U,y)
    # Second iteration
    r = d - np.dot(A,x)
    dz = solve(L,r)
    ddx = solve(U,dz)
    # Solution
    totx = x+ddx
    self.state.xhat_new = self.state.xhat+totx
    self.state.d2 = np.dot(totx.T,d)
    
  def invert(self,profiling=None):

    if ( profiling is not None ):
      if ( profiling == True ):
        pr = cProfile.Profile()
        pr.enable()

    # Assing values for the state vector
    xhat = self.apriori.x0

    # Iteration of the Newton-Gauss method to find the zero of the first
    #   derivative of the Gaussian PDF
    Iteration = 0

    # Initialize state vector per previous iteration to apriori.xa
    #   (same as apriori.X0)
    xhat_pre = self.apriori.xa

    fm = mirto_code_compute_F.radiance(self.cx)

    jj = self.cx.state_var_indx
    self.state.SaInv_ret = self.apriori.SaInv[jj,:]
    self.state.SaInv_ret = self.state.SaInv_ret[:,jj]
    self.state.xa = self.apriori.xa[jj]

    while (Iteration < self.cx.Iteration_limit):
      #
      # Compute F (forward model)
      #
      # print('... Compute_F')

      fm.compute_forward(xhat)
      fm.estimate_K(xhat)

      fm.compute_residuals(self.residuals)

      # Subselect from forward model output the channels used for the
      # inversion

      ii = self.cx.indx
      fm.K = fm.K[ii,:]
      fm.F = fm.F[ii]
      fm.wnF = fm.wnF[ii]
      fm.wnK = fm.wnK[ii]

      # Subselect from forward model output (Jacobians) and from
      # whole apriori the variables actually retrieved.
      fm.K = fm.K[:,jj]
      self.state.xhat = xhat[jj]
      self.state.xhat_pre = xhat_pre[jj]

      self.compute_chi_square()
      self.update_solution(fm)

      if ( Iteration > 0 ):
        ref_norm = min(self.hist.measurement_space_only)
        print ('Actual value of Norm   : ',self.norm.measurement_space_only)
        print ('Last low value of Norm : ',ref_norm)
        if (self.norm.measurement_space_only <= ref_norm):
          self.cx.gamma = self.cx.gamma/2.0
          xxdel = (100.0 * (ref_norm - self.norm.measurement_space_only) /
                            self.norm.measurement_space_only)
          print('UWPHYSRET residuals decreased by ',xxdel,'%')
        else:
          self.cx.gamma = self.cx.gamma*5.0
          xxdel = (100.0 * (self.norm.measurement_space_only - ref_norm) /
                            self.norm.measurement_space_only)
          print('UWPHYSRET residuals increased by ',xxdel,'%')

      xhat[jj] = self.state.xhat

      if (abs(self.state.d2) < self.converge):
        print('****  CONVERGED!!! ****')
        break
      else:
        if (Iteration < self.cx.Iteration_limit):
          # assign new value to the solution and continue the iterative solution
          # Note: Don't update xhat if this is the final iteration or it
          # will overwrite the solution
          xhat_pre[jj] = self.state.xhat
          xhat[jj] = self.state.xhat_new

      self.hist.measurement_space_only.append(self.norm.measurement_space_only)
      Iteration = Iteration + 1
      print ('Iteration = ', Iteration)
      print ('New Gamma = ', self.cx.gamma)
      print ('Distance  = ', abs(self.state.d2))
      print ('Wanted    = ', self.converge)

    if ( profiling is not None ):
      if ( profiling == True ):
        pr.disable()
        s = ()
        try:
          # Default to Python 3
          import io
          s = io.StringIO()
          ps = pstats.Stats(pr, stream=s)
          ps.strip_dirs().sort_stats('cumulative').print_stats()
          print(s.getvalue())
        except:
          # May be this is Python 2?
          import StringIO
          s = StringIO.StringIO()
          ps = pstats.Stats(pr, stream=s)
          ps.strip_dirs().sort_stats('cumulative').print_stats()
          print(s.getvalue())

    return(self.state)

if ( __name__ == '__main__' ):
  sys.path.append('/home/graziano/Software/pythoncode/oss')
  from oss4SHIS import oss4SHIS
  from os import path
  import time
  #
  # OSS init input
  #
  solar = 'solar_irradiances.nc'
  precomputed = 'leo.cris.0.05.nc'
  datapath = '/home/graziano/Software/pythoncode/data'
  oss = oss4SHIS(path.join(datapath,solar),path.join(datapath,precomputed))

  # This part of code must be repeated for each input profile in data
  # directory. Must find a way to have names here. Probably the errors
  # also can be preloaded.
  start = time.clock()
  profiling = False
  check_output = False
  inverter = mirto('/home/graziano/Software/pythoncode/data',oss)
  solution = inverter.invert(profiling)
  print('Elapsed Time in the Inversion: ',time.clock() - start,' s')

  if ( check_output ):
    print('Profiles')
    print('Pressure         Temperature             Water Vapor         O3')
    for i in range(0,61):
      print(inverter.cx.pressure_grid[i],solution.xhat[i],
            solution.xhat[i+61],solution.xhat[i+122])
    print('Skin temperature : ',solution.xhat[183])
    print('Surface Emissivity : ')
    for i in range(184,189):
      print(solution.xhat[i])
    print('Value of distance : ',solution.d2)
    try:
      import pylab as p
    except:
      print('Cannot use pylab. Is it installed?')
      sys.exit()
    x = solution.xhat[0:61]
    y = np.log(inverter.cx.pressure_grid[0:61])
    p.ylabel("Log Pressure")
    p.xlabel("Temperature")
    p.plot(x,y)
    p.plt.gca().invert_yaxis()
    p.show()
    p.ylabel("Pressure")
    p.xlabel("Mixing Ratio")
    x = np.exp(solution.xhat[61:122])
    y = inverter.cx.pressure_grid[0:61]
    p.plot(x,y)
    p.plt.gca().invert_yaxis()
    p.show()

