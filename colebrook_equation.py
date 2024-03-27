import numpy as np

# Inputs
e = 0.045E-3 # m
D = 0.289 # m
Re = 439280

def colebrook_equation(e, D, Re):

  # Function to calculate RHS of Colebrook Equation
  def RHS(e, D, Re, f):
      return -2.0*np.log10( (e/D)/3.7 + 2.51/(Re*np.sqrt(f)) )

  # Initial guess for f
  f = 0.01
  for iter in range(10):
      print(f'Iteration {iter}')
      colebrook_RHS = RHS(e, D, Re, f)
      f = 1/colebrook_RHS**2
      print(f'f = {f}')
  
  return f