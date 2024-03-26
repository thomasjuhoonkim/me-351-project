import numpy as np

# Inputs (from colebrook)
# D = 0.289 # m
# Re = 439280

# CONSTANTS
e = 0.0015E-3 # surface roughness (m)
V_rate = 1234567 # volumetric flow rate (m^3/s)

PIPE_TYPES = {
  'NPS 6': {
    'OD': 168.28,
    'WT': 7.112
  },
  'NPS 7': {
    'OD': 193.68,
    'WT': 7.645
  },
  'NPS 8': {
    'OD': 219.08,
    'WT': 8.179
  },
  'NPS 9': {
    'OD': 244.48,
    'WT': 8.687
  },
  'NPS 10': {
    'OD': 273.05,
    'WT': 9.271,
  },
  'NPS 12': {
    'OD': 323.85,
    'WT': 9.525,
  },
  'NPS 14': {
    'OD': 355.60,
    'WT': 9.525,
  },
  'NPS 16': {
    'OD': 406.40,
    'WT': 9.525,
  },
  'NPS 18': {
    'OD': 457.20,
    'WT': 9.525,
  },
  'NPS 20': {
    'OD': 508.00,
    'WT': 9.525,
  },
  'NPS 22': {
    'OD': 558.80,
    'WT': 9.525,
  },
  'NPS 24': {
    'OD': 609.60,
    'WT': 9.525,
  },
}

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

def get_inner_diameter(pipe_type):
  assert pipe_type in PIPE_TYPES, "Not a valid pipe type"

  OD, WT = 'OD', 'WT'
  return PIPE_TYPES[pipe_type][OD] - 2 * PIPE_TYPES[pipe_type][WT]

def get_pipe_cost(pipe_type, pipe_length):
  assert pipe_type in PIPE_TYPES, "Not a valid pipe type"

  OD = 'OD' # Nominal Diameter = Outer Diameter
  cost_per_length = 200 * (PIPE_TYPES[pipe_type][OD] / 150) ** 1.2
  return pipe_length * cost_per_length

def get_railway_crossing(x1, x1_pipe_type, x2):
  # railway crossing calculation

  # return cost, x2_pipe_type
  pass