#!/usr/bin/env python
# mirto_code_main.py
import numpy as np
from scipy.linalg import lu , solve
import cProfile, pstats, io
import mirto_code_configuration
import mirto_code_compute_F
import sys

# Empty class to store result
class mirto_state:
  pass

class mirto_norm:
  pass

class mirto:
  def __init__(self,datapath,oss):
    # Load configuration parameters and Input data for retrieval
    self.control = mirto_code_configuration.control(datapath,oss)
    self.obsErr = mirto_code_configuration.obsErr(self.control)
    self.apriori = mirto_code_configuration.apriori(self.control)
    self.state = mirto_state()

  def compute_chi_square(self,residuals):
    norm = mirto_norm()
    norm.measurement_space_only = (
      np.dot(
             np.dot(residuals.yobs_minus_yhat.T,self.obsErr.SeInv),
             residuals.yobs_minus_yhat/len(residuals.yobs_minus_yhat)) )
    norm.solution_space_only = (
      np.dot(
             np.dot(self.state.xhat-self.state.xhat_pre.T,
                    self.apriori.SaInv_ret),
             self.state.xhat-self.state.xhat_pre/len(self.state.xhat) ) )
    return(norm)

  def update_solution(self,fm,residuals):
    self.control.KtSeInv = np.dot(fm.K.T,self.obsErr.SeInv)
    self.control.KtSeInvK = np.dot(self.control.KtSeInv,fm.K)
    A = (self.control.KtSeInvK+(1.0+self.control.gamma)*self.apriori.SaInv_ret)
    dx = (self.state.xhat-self.state.xa)
    d = (np.dot(self.control.KtSeInv,residuals.yobs_minus_yhat) - 
         np.dot(self.apriori.SaInv_ret,dx))
    # Use iterative LU decomposition to determine the solution
    # First iteration
    (P,L,U) = lu(A)
    y = solve(L,d)
    x = solve(U,y)
    # Second iteration
    r = d - np.dot(A,x)
    dz = solve(L,r)
    ddx = solve(U,dz)
    # Solution
    self.state.xhat_new = x+ddx+self.state.xhat
    self.state.d2 = np.dot(self.state.xhat_new-self.state.xhat.T,d)
    
  def invert(self,profiling=None):
    if ( profiling is not None ):
      if ( profiling == True ):
        pr = cProfile.Profile()
        pr.enable()

    # Assing values for the state vector
    xhat = self.apriori.x0

    # Iteration of the Newton-Gauss method to find the zero of the first
    #   derivative of the Gaussian PDF
    Iteration = 1

    # Set initial values (very large) for d2
    self.control.d2 = 10.0**6

    # Initialize state vector per previous iteration to apriori.xa
    #   (same as apriori.X0)
    xhat_pre = self.apriori.xa

    fm = mirto_code_compute_F.radiance(self.control)

    while (Iteration <= self.control.Iteration_limit):
      #
      # Compute F (forward model)
      #
      # print('... Compute_F')
      fm.compute_forward(xhat)
      fm.estimate_K(xhat)
      residuals = fm.compute_residuals( )

      # Subselect from forward model output the channels used for the
      # inversion
      fm.K = fm.K[:,self.control.state_var_indx]
      fm.F = fm.F[self.control.indx]
      fm.wnF = fm.wnF[self.control.indx]
      fm.wnK = fm.wnK[self.control.indx]

      # Subselect from forward model output (Jacobians) and from
      # whole apriori the variables actually retrieved.
      fm.K = fm.K[self.control.indx,:]
      self.apriori.Sa_ret = self.apriori.Sa[self.control.state_var_indx,
                                            self.control.state_var_indx]
      self.apriori.SaInv_ret = self.apriori.SaInv[self.control.state_var_indx,
                                                  self.control.state_var_indx]
      self.state.xhat = xhat[self.control.state_var_indx]
      self.state.xhat_pre = xhat_pre[self.control.state_var_indx]
      self.state.xa = self.apriori.xa[self.control.state_var_indx]

      norm = self.compute_chi_square(residuals)

      self.update_solution(fm,residuals)

      sys.exit()

# UP TO HERE

      [state, self.control] = mirto_code_ml_update_xhat(state, fm, self.apriori,
                     self.obsErr, residuals, self.control)

      historical.measurement_space_only[Iteration] = norm.measurement_space_only
      historical.d2[Iteration] = state.d2

      if (Iteration > 1):
        ref_norm = min(historical.measurement_space_only[0:end-1]);
        if  (norm.measurement_space_only <= ref_norm):
          print('Residual decreasing, solution saved')
          self.control.gamma = self.control.gamma/2
          xxdel = (100.0 * (ref_norm - norm.measurement_space_only) /
                            norm.measurement_space_only)
          print('UWPHYSRET residuals decreased by ',xxdel)
        else:
          self.control.gamma = self.control.gamma*5.0
          xxdel = (100.0 * (norm.measurement_space_only - ref_norm) /
                            norm.measurement_space_only)
          print('UWPHYSRET residuals increased by ',xxdel)
          print('Residual increasing, solution NOT saved')
          his_d2 = min(abs(historical.d2[0:end-1]))
      else:
        ref_d2 = 10.0**10
        ref_norm = 10.0**15

      xhat[self.control.state_var_indx] = state.xhat

      if (abs(state.d2) < 0.03 * length(self.control.state_var_indx)):
        print('****  CONVERGED!!! ****')
        break
      else:
        print('****  Convergence test failed.  ')
        if (Iteration < self.control.Iteration_limit):
          # assign new value to the solution and continue the iterative solution
          # Note: Don't update xhat if this is the final iteration or it
          # will overwrite the solution
          print('... Update_xhat')
          xhat_pre[self.control.state_var_indx] = state.xhat
          xhat[self.control.state_var_indx] = state.xhat_new
      Iteration = Iteration+1

    return(state)

    if ( profiling is not None ):
      if ( profiling == True ):
        pr.disable()
        s = io.StringIO()
        sortby = 'cumulative'
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        print(s.getvalue())

if ( __name__ == '__main__' ):
  sys.path.append('/home/graziano/Software/pythoncode/oss')
  from oss4SHIS import oss4SHIS
  from os import path
  #
  # OSS init input
  #
  solar = 'solar_irradiances.nc'
  precomputed = 'leo.cris.0.05.nc'
  datapath = '/home/graziano/Software/pythoncode/data'
  oss = oss4SHIS(path.join(datapath,solar),path.join(datapath,precomputed))
  inverter = mirto('/home/graziano/Software/pythoncode/data',oss)
  solution = inverter.invert()
