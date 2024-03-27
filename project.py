import numpy as np
import csv
import collections
import bisect

# ========== CONSTANTS START ==========

V_RATE = 0.05 # volumetric flow rate (m^3/s)

E = 0.0015e-3 # surface roughness (m)
P = 1000 # density (kg/m^2)
MU = 8.9e-4 # dynamic viscosity (Pa.s)
G = 9.81 # gravity (N/kg)

TOTAL_DISTANCE = 62550 # distance from lake ontario to E5 (m)

# pressure constants
ATM_PRESSURE = 101325 # atmospheric pressure (Pa)
MIN_PRESSURE = 275000 + ATM_PRESSURE # 275 kPag is gauge pressure (Pa)
MAX_PRESSURE = 700000 + ATM_PRESSURE # 700 kPag is gauge pressure (Pa)

# minor loss coefficients (dimensionless)
KL_180_DEG_RETURN_BEND = 0.2
KL_90_DEG_MITRE_BEND_VANES = 0.2
KL_90_DEG_SMOOTH_BEND = 0.3
KL_GATE_VALVE = 0.2
KL_FLOW_METER = 7
KL_SWING_CHECK_VALVE = 2

# NPS pipe types and their properties
#   OD = Outer Diameter (mm)
#   WT = Wall Thickness (mm)
#   DN = Nominal DN (dimensionless)
PIPE_TYPES = {
  'NPS 6': {
    'OD': 168.28,
    'WT': 7.112,
    'DN': 150,
  },
  'NPS 8': {
    'OD': 219.08,
    'WT': 8.179,
    'DN': 200,
  },
  'NPS 10': {
    'OD': 273.05,
    'WT': 9.271,
    'DN': 250,
  },
  'NPS 12': {
    'OD': 323.85,
    'WT': 9.525,
    'DN': 300,
  },
  'NPS 14': {
    'OD': 355.60,
    'WT': 9.525,
    'DN': 350,
  },
  'NPS 16': {
    'OD': 406.40,
    'WT': 9.525,
    'DN': 400,
  },
  'NPS 18': {
    'OD': 457.20,
    'WT': 9.525,
    'DN': 450,
  },
  'NPS 20': {
    'OD': 508.00,
    'WT': 9.525,
    'DN': 500,
  },
  'NPS 22': {
    'OD': 558.80,
    'WT': 9.525,
    'DN': 550,
  },
  'NPS 24': {
    'OD': 609.60,
    'WT': 9.525,
    'DN': 600,
  },
}

# costs of common components of the pipe system
COMPONENT_COSTS = {
  'VALVE': 50000, # per valve
  'PUMP STATION': 1500000, # per station (including valves)
  'RAILWAY CROSSING': 100000, # per crossing
  'ROAD CROSSING': 175000, # per crossing
  'WATER BODY CROSSING': 300000, # per crossing
}

# ========== CONSTANTS END ==========



# ========== VARIABLE CALCULATION FUNCTIONS START ==========

def get_friction_factor_colebrook(e, D, Re):
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

def get_reynolds_number(p, v_avg, D, mu):
  return p * v_avg * D / mu

def get_inner_diameter(pipe_type):
  assert pipe_type in PIPE_TYPES, "Not a valid pipe type"
  inner_diameter_mm = PIPE_TYPES[pipe_type]['OD'] - 2 * PIPE_TYPES[pipe_type]['WT']
  return inner_diameter_mm * 1e-3 # convert to meters

def get_cross_sectional_area(pipe_type):
  assert pipe_type in PIPE_TYPES, "Not a valid pipe type"
  return np.pi / 4 * get_inner_diameter(pipe_type) ** 2

def get_average_velocity(pipe_type):
  assert pipe_type in PIPE_TYPES, "Not a valid pipe type"
  return V_RATE / get_cross_sectional_area(pipe_type)

def get_pipe_cost(pipe_type, pipe_length):
  assert pipe_type in PIPE_TYPES, "Not a valid pipe type"
  cost_per_meter = 200 * ((PIPE_TYPES[pipe_type]['DN'] / 150) ** 1.2)
  return pipe_length * cost_per_meter

# ========== VARIABLE CALCULATION FUNCTIONS END ==========



# ========== PIPE SECTIONS START ==========

def add_regular_pipe_section(x1, y1, x2, y2):
  pass

def get_railway_crossing(x1, pipe_type_1, y1, x2, x2_pipe_type, y2):
  pass

# pump station can be treated as having no length
# pump station has no elevation changes (zero hydrostatic pressure changes)
# pump station has no internal major losses
# pump station is a single minor loss
# pump station supplies max allowable pressure minus minor losses
def add_pump_station(x, pipe_type, valve_locations):

  NUM_90_DEG_SMOOTH_BENDS = 4
  NUM_GATE_VALVES = 6
  NUM_FLOW_METERS = 1
  NUM_SWING_CHECK_VALVES = 1

  # check for redundant valves (more than 2 valves every 300 m) and if present, remove them
  while len(valve_locations) >= 2 and x - valve_locations[-2] <= 300:
    valve_locations.pop()

  for _ in range(NUM_GATE_VALVES):
    valve_locations.append(x)

  KL_COMBINED = NUM_90_DEG_SMOOTH_BENDS * KL_90_DEG_SMOOTH_BEND + \
                NUM_GATE_VALVES * KL_GATE_VALVE + \
                NUM_FLOW_METERS * KL_FLOW_METER + \
                NUM_SWING_CHECK_VALVES * KL_SWING_CHECK_VALVE

  minor_losses = KL_COMBINED * P * (get_average_velocity(pipe_type) ** 2) / 2

  return MAX_PRESSURE - minor_losses

# ========== PIPE SECTIONS END ==========



# ========== UTILITY FUNCTIONS START ==========

def remove_within_range_list_of_tuples(l, lower, upper):
  return [tup for tup in l if tup[0] <= lower or upper <= tup[0]]

def custom_round(x, base=50):
  return base * round(x/base)

# import the csv into the script ('name' is for logging purposes)
def import_data(name, path):
  print(f'Importing {name}...')
  data = collections.deque()

  with open(path) as csv_file:
    csv_reader = csv.reader(csv_file)
    line_count = 0
    for row in csv_reader:
      if line_count == 0:
        print(f'Column names are "{", ".join(row)}"')
        line_count += 1
      else:
        row = tuple(float(num) for num in row)
        data.append(row)
        line_count += 1
    print(f'Processed {line_count} lines')

  print("-----")
  return data

def export_data(name, path):
  pass

# ========== UTILITY FUNCTIONS END ==========



# ========== SCRIPT START ==========

if __name__ == "__main__":
  # import data
  print("===== IMPORT DATA =====")
  ROAD_CROSSINGS = import_data('Road Crossings', 'data/road_crossings.csv')
  RAILWAY_CROSSINGS = import_data('Railway Crossings', 'data/railway_crossings.csv')
  WATER_CROSSINGS = import_data('Water Crossings', 'data/water_crossings.csv')
  ELEVATION_PROFILE = import_data('Elevation Profile', 'data/elevation_profile.csv')


  # initialize runtime variables
  valve_locations = [] # array of all valve locations x-coor
  crossing_valve_locations = [] # array of pairs of valves at start and end of crossings
  pressure = 0
  cost = 0
  num_pump_stations = 0


  # determine crossing valve locations and add crossing elevations
  print("===== DETERMINE CROSSING VALVE LOCATIONS AND ADD CROSSING LEVATIONS =====")
  for x, y in ROAD_CROSSINGS:
    x, y = float(x), float(y)

    # elevation of the crossing might not be relevant
    bisect.insort_left(ELEVATION_PROFILE, (x, y), key=lambda x: x[0])

  for x, y in RAILWAY_CROSSINGS:
    x, y = float(x), float(y)

    x_entry = max(0, x - 100) # ensure location is not less than 0
    x_exit = min(TOTAL_DISTANCE, x + 100) # ensure exit valve is not placed after E5

    bisect.insort_left(crossing_valve_locations, (x_entry, x_exit, 'RAILWAY CROSSING'), key=lambda x: x[0])

    # elevation of the crossing might not be relevant here (remove this for cost savings)
    bisect.insort_left(ELEVATION_PROFILE, (x, y), key=lambda x: x[0])

  for x1, y1, x2, y2 in WATER_CROSSINGS:
    x1, y1, x2, y2 = float(x1), float(y1), float(x2), float(y2)

    x_entry = max(0, x1 - 100) # ensure location is not less than 0
    x_exit = min(TOTAL_DISTANCE, x2 + 100) # ensure exit valve is not placed after E5

    bisect.insort_left(crossing_valve_locations, (x_entry, x_exit, 'WATER CROSSING'), key=lambda x: x[0])

    bisect.insort_left(ELEVATION_PROFILE, (x1, y1), key=lambda x: x[0])
    bisect.insort_left(ELEVATION_PROFILE, (x2, y2), key=lambda x: x[0])

    # get rid of elevation profiles within water bodies
    ELEVATION_PROFILE = remove_within_range_list_of_tuples(ELEVATION_PROFILE, x1, x2)
  

  # determine all valve locations
  print("===== DETERMINE ALL VALVE LOCATIONS =====")
  # add placeholder elements
  crossing_valve_locations.insert(0, (0, 0, 'START'))
  crossing_valve_locations.append((TOTAL_DISTANCE, TOTAL_DISTANCE, 'END'))

  # sliding window algorithm
  k = 2 # window size
  for i in range(len(crossing_valve_locations) - k + 1):
    crossing1 = crossing_valve_locations[i]
    crossing2 = crossing_valve_locations[i+k-1]

    distance = crossing2[0] - crossing1[1]
    num_sections = int(np.ceil(distance / 300))
    print(crossing1, crossing2)
    print(num_sections, distance)

    # number of separations is equal to num_valves_between - 1
    num_valves_between = num_sections - 1
    for n in range(1, num_valves_between + 1):
      distance_between_valves = distance / num_sections
      valve_locations.append(crossing1[0] + n * distance_between_valves)
    
    # append the valves at the crossings excluding the first and last placeholder elements
    if i < len(crossing_valve_locations) - k + 1:
      valve_locations.append(crossing2[0])
      valve_locations.append(crossing2[1])
  
  



  # ---------- ALL POSSIBLE REYNOLDS NUMBERS BY VARYING PIPE SIZE ----------
  # this shows that no matter the pipe size, the flow will
  # always remain turbulent as the reynold's number will
  # ALWAYS be greater than 4000

  # for pipe_type in PIPE_TYPES.keys():
  #   ID = get_inner_diameter(pipe_type)
  #   CA = get_cross_sectional_area(pipe_type)
  #   Re = get_reynolds_number(P, V_RATE / CA, ID, MU)
  #   print(Re)
   

# ========== SCRIPT END ==========