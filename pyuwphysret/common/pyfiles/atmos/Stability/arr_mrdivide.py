#!/usr/bin/env pyton
import numpy as np

def arr_arr_mrdivide(a,b):
  M = np.vstack([a, np.zeros(len(a))]).T
  x = np.linalg.lstsq(M, b)[0][0]
  return x

def arr_val_mrdivide(a,b):
  B = np.zeros(len(a))+b
  M = np.vstack([a, np.zeros(len(a))]).T
  x = np.linalg.lstsq(M, B)[0][0]
  return x

def val_arr_mrdivide(a,b):
  A = np.zeros(len(b))+a
  M = np.vstack([A, np.zeros(len(A))]).T
  x = np.linalg.lstsq(M, b)[0][0]
  return x

if __name__ == '__main__':
  a = np.array(( 5.,  6.,  7.,  8. ))
  b = np.array(( 1.,  1.,  1.,  1. ))
  print(arr_arr_mrdivide(a,b))
  print(arr_arr_mrdivide(b,a))
  c = 1.0
  print(arr_val_mrdivide(a,c))
  b = a
  a = 1.0
  print(val_arr_mrdivide(a,b))
