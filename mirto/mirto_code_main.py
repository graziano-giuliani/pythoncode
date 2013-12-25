#!/usr/bin/env python
# mirto_code_main.py
import numpy as num
import cProfile, pstats, io
import mirto_code_configuration
import sys

def mirto_code_main(datapath,oss,profiling=None):
  """solution = mirto_code_main

This function is part of the Matlab implementation of the

            MTG-IRS Real Time Operational code

The function is the main driver for the inversion process

Input: None

Output: Completion code

Calling sequence:

Developed by: P. Anontelli, G. Giuliani, T. Cherubini

Thu 12 dec 2013, 10.32.53, CET
  """

  if ( profiling is not None ):
    if ( profiling == True ):
      pr = cProfile.Profile()
      pr.enable()

  # Load configuration parameters and Input data for retrieval
  control = mirto_code_configuration.control(datapath,oss)
  obsErr = mirto_code_configuration.obsErr(control)
  apriori = mirto_code_configuration.apriori(control)

  # Assing values for the state vector
  xhat = apriori.x0

  # Iteration of the Newton-Gauss method to find the zero of the first
  #   derivative of the Gaussian PDF
  Iteration = 1

  # Set initial values (very large) for d2
  control.d2 = 10.0**6

  # Initialize state vector per previous iteration to apriori.xa
  #   (same as apriori.X0)
  xhat_pre = apriori.xa

  fm = radiance(control)

  while (Iteration <= control.Iteration_limit):
    #
    #
    # Compute F (forward model)
    #
    #
    # print('... Compute_F')
    fm.compute(xhat)
    fm.estimate_K(xhat)

# UP TO HERE

    residuals = mirto_code_compute_residuals(control,fm)

    fm = mirto_code_extract_selected_channels(control,fm)

    [apriori, fm, state] = mirto_code_extract_selected_variables(apriori,
                     control, fm, xhat, xhat_pre)

    norm = mirto_code_compute_chi_square(state, residuals, obsErr, apriori)

    [state, control] = mirto_code_ml_update_xhat(state, fm, apriori,
                   obsErr, residuals, control)

    historical.measurement_space_only[Iteration] = norm.measurement_space_only
    historical.d2[Iteration] = state.d2

    if (Iteration > 1):
      ref_norm = min(historical.measurement_space_only[0:end-1]);
      if  (norm.measurement_space_only <= ref_norm):
        print('Residual decreasing, solution saved')
        control.gamma = control.gamma/2
        xxdel = (100.0 * (ref_norm - norm.measurement_space_only) /
                          norm.measurement_space_only)
        print('UWPHYSRET residuals decreased by ',xxdel)
      else:
        control.gamma = control.gamma*5.0
        xxdel = (100.0 * (norm.measurement_space_only - ref_norm) /
                          norm.measurement_space_only)
        print('UWPHYSRET residuals increased by ',xxdel)
        print('Residual increasing, solution NOT saved')
        his_d2 = min(abs(historical.d2[0:end-1]))
    else:
      ref_d2 = 10.0**10
      ref_norm = 10.0**15

    xhat[control.state_var_indx] = state.xhat

    if (abs(state.d2) < .03 * length(control.state_var_indx)):
      print('****  CONVERGED!!! ****')
      break
    else:
      print('****  Convergence test failed.  ')
      if (Iteration < control.Iteration_limit):
        # assign new value to the solution and continue the iterative solution
        # Note: Don't update xhat if this is the final iteration or it
        # will overwrite the solution
        print('... Update_xhat')
        xhat_pre[control.state_var_indx] = state.xhat
        xhat[control.state_var_indx] = state.xhat_new

    Iteration = Iteration+1
    return(state)

  if ( profiling is not None ):
    if ( profiling == True ):
      pr.disable()
      s = io.StringIO.StringIO()
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
  solution = mirto_code_main('/home/graziano/Software/pythoncode/data',oss)
