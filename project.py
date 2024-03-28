import numpy as np
import csv
import collections
import bisect

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

def export_data(name, path, data):
  print(f'Exporting {name}...')

  with open(path, mode='w') as csv_file:
    csv_writer = csv.writer(csv_file)
    line_count = 0
    for row in data:
      csv_writer.writerow(row)
      line_count += 1
    print(f'Processed {line_count} lines')
  
  print("-----")
  return None

# ========== UTILITY FUNCTIONS END ==========



# ========== CONSTANTS START ==========

V_RATE = 0.05 # volumetric flow rate (m^3/s)

E = 0.0015e-3 # surface roughness (m)
P = 1000 # density (kg/m^2)
MU = 8.9e-4 # dynamic viscosity (Pa.s)
G = 9.81 # gravity (N/kg)
ALPHA = 1.05 # kinetic energy correction factor for fully turbulent flow

MAX_DISTANCE_BETWEEN_VALVES = 300 # distance between valves (m)
TOTAL_DISTANCE = 62550 # distance from lake ontario to E5 (m)
LAKE_ONTARIO_PIPE_START_DEPTH = 40 # depth of pipe in lake ontario at the start location (m)

# pressure constants
PRESSURE_ATM = 101325 # atmospheric pressure (Pa)
PRESSURE_MIN = 275000 + PRESSURE_ATM # 275 kPag is in gauge pressure (Pa)
PRESSURE_MAX = 700000 + PRESSURE_ATM # 700 kPag is in gauge pressure (Pa)

# minor loss coefficients (dimensionless)
KL_PIPE_OUTLET = ALPHA
KL_PIPE_INLET_PROTRUDED = 0.8
KL_180_DEG_RETURN_BEND = 0.2
KL_90_DEG_MITRE_BEND_VANES = 0.2
KL_90_DEG_SMOOTH_BEND = 0.3
KL_GATE_VALVE = 0.2
KL_FLOW_METER = 7
KL_SWING_CHECK_VALVE = 2

# component types
ROAD_CROSSING = "ROAD CROSSING"
RAILWAY_CROSSING = "RAILWAY CROSSING"
WATER_CROSSING = "WATER CROSSING"
VALVE = "VALVE"
PUMP_STATION = "PUMP STATION"
STANDARD_PIPE_SECTION = "STANDARD_PIPE_SECTION"

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
  VALVE: 50000, # per valve
  PUMP_STATION: 1500000, # per station (including valves)
  RAILWAY_CROSSING: 100000, # per crossing
  ROAD_CROSSING: 175000, # per crossing
  WATER_CROSSING: 300000, # per crossing
}

# ========== CONSTANTS END ==========



# ========== VARIABLE CALCULATION FUNCTIONS START ==========

def get_inner_diameter(pipe_type):
  assert pipe_type in PIPE_TYPES, "Not a valid pipe type"
  inner_diameter_mm = PIPE_TYPES[pipe_type]['OD'] - 2 * PIPE_TYPES[pipe_type]['WT']
  return inner_diameter_mm * 1e-3 # convert to meters

def get_cross_sectional_area(inner_diameter):
  return np.pi / 4 * inner_diameter ** 2

def get_average_velocity(V_rate, area):
  return V_rate / area

def get_pipe_cost(pipe_type, pipe_length):
  assert pipe_type in PIPE_TYPES, "Not a valid pipe type"
  cost_per_meter = 200 * ((PIPE_TYPES[pipe_type]['DN'] / 150) ** 1.2)
  return pipe_length * cost_per_meter

def get_reynolds_number(p, v_avg, D, mu):
  return p * v_avg * D / mu

# directly from the project materials
def get_friction_factor_colebrook(e, D, Re):
  # Function to calculate RHS of Colebrook Equation
  def RHS(e, D, Re, f):
      return -2.0*np.log10( (e/D)/3.7 + 2.51/(Re*np.sqrt(f)) )

  # Initial guess for f
  f = 0.01
  for iter in range(10):
      # print(f'Iteration {iter}')
      colebrook_RHS = RHS(e, D, Re, f)
      f = 1/colebrook_RHS**2
      # print(f'f = {f}')
  
  return f

def get_major_loss(f, L, D, p, v_avg):
  return f * (L / D) * (p * (v_avg ** 2) / 2)

def get_minor_loss(KL, p, v_avg):
  return KL * (p * (v_avg ** 2) / 2)

def get_hydrostatic_pressure(p, g, z):
  return p * g * z

def get_kilo_pascals(pressure):
  return pressure / 1000

# ========== VARIABLE CALCULATION FUNCTIONS END ==========



# ========== PIPE SECTIONS START ==========

# 1 ╗
#   ║ ↑           Major losses along length of pipes
#   ║ y1-y2       Minor losses in two bends
#   ║ ↓           2 x 90 deg mitre bends with vanes
#   ║             KL_EQUIVALENT= 2 * 0.2
#   ║ ←x1-x2→     P2 - P1 = pg(z1-z2) - dP_losses
#   ╚═════════ 2
# NOTE: This is "dy" sensitive meaning the section will flip
#       vertically if 1 is lower than 2 (climbing pipe section)
# NOTE: This is an approximation of upward and downward
#       elevation changes, in real life, if the pipes are
#       implemented as such, some sections will protrude
# NOTE: This approximation can be improved using a different
#       approximation function for the pipe section
#       make sure to modify the function being used in the
#       script as defined in the variable "add_pipe_section"
def add_90_deg_pipe_section(pipe_type, dx, dy):

  ID = get_inner_diameter(pipe_type)
  A = get_cross_sectional_area(ID)
  V_AVG = get_average_velocity(V_RATE, A)
  RE = get_reynolds_number(P, V_AVG, ID, MU)
  F = get_friction_factor_colebrook(E, ID, RE)

  NUM_90_DEG_MITRE_BENDS_VANES = 2

  KL_EQUIVALENT = NUM_90_DEG_MITRE_BENDS_VANES * KL_90_DEG_MITRE_BEND_VANES

  # dx = x1 - x2
  # dy = y1 - y2
  L_vertical = abs(dy)
  L_horizontal = abs(dx)

  # major losses
  dP_major_vertical = get_major_loss(F, L_vertical, ID, P, V_AVG)
  dP_major_horizontal = get_major_loss(F, L_horizontal, ID, P, V_AVG)

  # minor losses
  dP_minor_bends = get_minor_loss(KL_EQUIVALENT, P, V_AVG)

  # hydrostatic pressure
  dP_hydrostatic = get_hydrostatic_pressure(P, G, dy)

  # total losses
  dP_losses_total = dP_major_vertical + dP_major_horizontal + dP_minor_bends

  # pressure drop
  dP = dP_hydrostatic - dP_losses_total

  # cost
  pipe_cost_vertical = get_pipe_cost(pipe_type, L_vertical)
  pipe_cost_horizontal = get_pipe_cost(pipe_type, L_horizontal)

  # total cost
  pipe_cost_total = pipe_cost_vertical + pipe_cost_horizontal

  return dP, pipe_cost_total

# using the add_section as a parameter fed function so that
# we can change the implementation of building a pipe
# between two points for potential cost savings
def add_road_crossing(dx, dy, pipe_type, add_pipe_section=add_90_deg_pipe_section):
  pressure_delta, cost = add_pipe_section(pipe_type, dx, dy)
  return pressure_delta, cost + COMPONENT_COSTS[ROAD_CROSSING]

# using the add_section as a parameter fed function so that
# we can change the implementation of building a pipe
# between two points for potential cost savings
def add_railway_crossing(dx, dy, pipe_type, add_pipe_section=add_90_deg_pipe_section):
  pressure_delta, cost = add_pipe_section(pipe_type, dx, dy)
  return pressure_delta, cost + COMPONENT_COSTS[RAILWAY_CROSSING]

# using the add_section as a parameter fed function so that
# we can change the implementation of building a pipe
# between two points for potential cost savings
def add_water_crossing(dx, dy, pipe_type, add_pipe_section=add_90_deg_pipe_section):
  pressure_delta, cost = add_pipe_section(pipe_type, dx, dy)
  return pressure_delta, cost + COMPONENT_COSTS[WATER_CROSSING]

# pump station can be treated as having no length
# pump station has no elevation changes (zero hydrostatic pressure changes)
# pump station has no internal major losses
# pump station is a single minor loss
# pump station supplies max allowable pressure minus minor losses
def add_pump_station(pipe_type):

  ID = get_inner_diameter(pipe_type)
  A = get_cross_sectional_area(ID)

  NUM_90_DEG_SMOOTH_BENDS = 4
  NUM_GATE_VALVES = 6
  NUM_FLOW_METERS = 1
  NUM_SWING_CHECK_VALVES = 1

  KL_EQUIVALENT = NUM_90_DEG_SMOOTH_BENDS * KL_90_DEG_SMOOTH_BEND + \
                NUM_GATE_VALVES * KL_GATE_VALVE + \
                NUM_FLOW_METERS * KL_FLOW_METER + \
                NUM_SWING_CHECK_VALVES * KL_SWING_CHECK_VALVE

  minor_losses = KL_EQUIVALENT * P * (get_average_velocity(V_RATE, A) ** 2) / 2

  return PRESSURE_MAX - minor_losses, COMPONENT_COSTS[PUMP_STATION]

def add_valve(pipe_type):

  ID = get_inner_diameter(pipe_type)
  A = get_cross_sectional_area(ID)
  V_AVG = get_average_velocity(V_RATE, A)

  # minor losses of valve
  dP_valve = -1 * get_minor_loss(KL_GATE_VALVE, P, V_AVG)

  return dP_valve, COMPONENT_COSTS[VALVE]

# ========== PIPE SECTIONS END ==========








# ========== SCRIPT START ==========

if __name__ == "__main__":
  # import data
  print() # spacer for legibility
  print("===== IMPORT DATA =====")
  ROAD_CROSSINGS = import_data('Road Crossings', 'input_data/road_crossings.csv')
  RAILWAY_CROSSINGS = import_data('Railway Crossings', 'input_data/railway_crossings.csv')
  WATER_CROSSINGS = import_data('Water Crossings', 'input_data/water_crossings.csv')
  ELEVATION_PROFILE = import_data('Elevation Profile', 'input_data/elevation_profile.csv')


  # initialize runtime variables - valve locations
  valve_locations = [] # array of all valve locations x-coor
  crossing_locations = []
  crossing_valve_locations = [] # array of pairs of valves at start and end of crossings


  # determine crossing valve locations
  print() # spacer for legibility
  print("===== DETERMINE VALVES ON EITHER END OF CROSSINGS =====")
  for x, y in ROAD_CROSSINGS:
    x, y = float(x), float(y)

    bisect.insort_left(crossing_locations, (x, x, ROAD_CROSSING), key=lambda x: x[0])


  for x, y in RAILWAY_CROSSINGS:
    x, y = float(x), float(y)

    bisect.insort_left(crossing_locations, (x, x, RAILWAY_CROSSING), key=lambda x: x[0])

    x_entry = max(0, x - 100) # ensure location is not less than 0
    x_exit = min(TOTAL_DISTANCE, x + 100) # ensure exit valve is not placed after E5

    bisect.insort_left(crossing_valve_locations, (x_entry, x_exit, RAILWAY_CROSSING), key=lambda x: x[0])


  for x1, y1, x2, y2 in WATER_CROSSINGS:
    x1, y1, x2, y2 = float(x1), float(y1), float(x2), float(y2)

    bisect.insort_left(crossing_locations, (x1, x2, WATER_CROSSING), key=lambda x: x[0])

    x_entry = max(0, x1 - 100) # ensure location is not less than 0
    x_exit = min(TOTAL_DISTANCE, x2 + 100) # ensure exit valve is not placed after E5

    bisect.insort_left(crossing_valve_locations, (x_entry, x_exit, WATER_CROSSING), key=lambda x: x[0])

    bisect.insort_left(ELEVATION_PROFILE, (x1, y1), key=lambda x: x[0])
    bisect.insort_left(ELEVATION_PROFILE, (x2, y2), key=lambda x: x[0])

    # get rid of elevation profiles within water bodies
    ELEVATION_PROFILE = remove_within_range_list_of_tuples(ELEVATION_PROFILE, x1, x2)

  for crossing_pairs in crossing_valve_locations:
    print(crossing_pairs)


  # determine all valve locations
  print() # spacer for legibility
  print("===== DETERMINE ALL VALVE LOCATIONS BETWEEN CROSSINGS =====")
  # add placeholder elements
  crossing_valve_locations.insert(0, (0, 0, 'START'))
  crossing_valve_locations.append((TOTAL_DISTANCE, TOTAL_DISTANCE, 'END'))


  # sliding window algorithm
  # there is no valve requirements for roads
  k = 2 # window size
  for i in range(len(crossing_valve_locations) - k + 1):
    crossing1 = crossing_valve_locations[i]
    crossing2 = crossing_valve_locations[i+k-1]

    distance = crossing2[0] - crossing1[1]
    num_sections = int(np.ceil(distance / MAX_DISTANCE_BETWEEN_VALVES))

    # number of sections is equal to num_valves_between + 1
    # example
    # ---- | ---- | ----
    # 3 sections, 2 valves
    num_valves_between = num_sections - 1
    for n in range(1, num_valves_between + 1):
      distance_between_valves = distance / num_sections
      valve_locations.append(crossing1[1] + n * distance_between_valves)
    
    # append the valves at the crossings excluding the first and last placeholder elements
    if i < len(crossing_valve_locations) - k:
      valve_locations.append(crossing2[0])
      valve_locations.append(crossing2[1])
    
    print(f"{num_valves_between} valves spaced {distance / num_sections} meters for {distance} meters between {crossing1[2]} at {crossing1[1]} and {crossing2[2]} at {crossing2[0]}")


  # initialize runtime variables - cost calculation
  minimum_cost = float('inf')
  minimum_cost_pipe_type = ''
  minimum_cost_pressure_data = [] # (x, pressure)
  minimum_cost_pump_station_locations = [] # (x, elevation)
  

  # determine the cheapest cost by varying pipe size
  print() # spacer for legibility
  print("===== VARY PIPE SIZE AND GET CHEAPEST COST (START AT LAKE ONTARIO) =====")
  for pipe_type in PIPE_TYPES.keys():
    
    print() # spacer for legibility
    print(f"----- Calculation for {pipe_type} -----")

    # initialize constants
    ID = get_inner_diameter(pipe_type)
    A = get_cross_sectional_area(ID)
    V_AVG = get_average_velocity(V_RATE, A)

    # initialize local variables
    pressure = get_hydrostatic_pressure(P, G, LAKE_ONTARIO_PIPE_START_DEPTH) + PRESSURE_ATM
    total_cost = 0
    crossing_index = 0
    valve_index = 0
    pressure_data = [] # (x, pressure)
    pump_station_locations = [] # (x, elevation)

    # initialize standard pipe function
    add_pipe_section = add_90_deg_pipe_section

    # based on the absolute pressure at the pipe inlet
    # when submerged under lake ontario, the pipe will
    # always begin with around 493 kPa
    print("PRESSURE AT INLET:", get_kilo_pascals(pressure), "kPa")

    # minor loss in pipe inlet at start
    pressure -= get_minor_loss(KL_PIPE_INLET_PROTRUDED, P, V_AVG)

    # pipe protrusion length cost
    PIPE_PROTRUSION_LENGTH = 0.1 * get_inner_diameter(pipe_type)
    total_cost += get_pipe_cost(pipe_type, PIPE_PROTRUSION_LENGTH)

    # sliding window algorithm - defines each individual section between x1 and x2
    k = 2 # window size
    for i in range(len(ELEVATION_PROFILE) - k + 1):

      x1, y1 = ELEVATION_PROFILE[i]
      x2, y2 = ELEVATION_PROFILE[i+k-1]

      dx = x1 - x2
      dy = y1 - y2

      pressure_delta = 0
      section_cost = 0
      is_crossing = False

      # log pressure data
      pressure_data.append((x1, pressure))

      # check if starting pressure is enough, if not enough we know
      # that pressure went below minimum pressure in the previous
      # section despite adding a pump, therefore, this pipe fails
      if pressure < PRESSURE_MIN:
        print(f"XXXXX {pipe_type} TERMINATED DUE TO INSUFFICIENT PRESSURE IN THE PREVIOUS SECTION DESPITE A PUMP ({get_kilo_pascals(pressure)} kPa < {get_kilo_pascals(PRESSURE_MIN)} kPa) XXXXX")
        break

      # check crossing index
      if crossing_index < len(crossing_locations):

        crossing_data = crossing_locations[crossing_index]
        crossing_x1 = crossing_data[0]
        crossing_x2 = crossing_data[1]
        crossing_type = crossing_data[2]

        # check if at crossing
        if x1 <= crossing_x1 <= x2 and x1 <= crossing_x2 <= x2:

          isCrossing = True
          
          if crossing_type == ROAD_CROSSING:
            print(f"Adding {crossing_type} at {crossing_x1} between {x1} and {x2}")
            pressure_delta, section_cost = add_road_crossing(dx, dy, pipe_type, add_pipe_section)
            
          if crossing_type == RAILWAY_CROSSING:
            print(f"Adding {crossing_type} at {crossing_x1} between {x1} and {x2}")
            pressure_delta, section_cost = add_railway_crossing(dx, dy, pipe_type, add_pipe_section)
            
          if crossing_type == WATER_CROSSING:
            print(f"Adding {crossing_type} at {crossing_x1} and {crossing_x2} between {x1} and {x2}")
            pressure_delta, section_cost = add_water_crossing(dx, dy, pipe_type, add_pipe_section)

          crossing_index += 1
      
      # standard pipe section if not a crossing
      if not is_crossing:
        pressure_delta, section_cost = add_pipe_section(pipe_type, dx, dy)

        # check if valves are within range, if so, add minor losses to pressure_delta
        check_index = valve_index
        num_valves_in_range = 0
        while check_index < len(valve_locations) and x1 <= valve_locations[check_index] <= x2:
          num_valves_in_range += 1
          check_index += 1
        
        # if there are valves within the range, add them to pressure drop and cost
        if num_valves_in_range > 0:
          for _ in range(num_valves_in_range):
            pressure_delta_valve, valve_cost = add_valve(pipe_type)
            pressure_delta += pressure_delta_valve
            section_cost += valve_cost

      # add a pump station at x1 if pressure is not sufficient
      if pressure + pressure_delta < PRESSURE_MIN:
        print(f"Adding {PUMP_STATION} at {x1} due to insufficient pressure ({get_kilo_pascals(pressure + pressure_delta)} kPa < {PRESSURE_MIN} kPa)")
        pressure, pump_cost = add_pump_station(pipe_type)
        section_cost += pump_cost

        # log pump station data
        pump_station_locations.append((x1, y1))

        # if a pump is added, we have to remove redundant valves
        # in both forward and backward directions to save money
        # we look at one valve beyond the nearest valve and see if
        # that distance is less than or equal to 300 meters
        # meaning that the valve nearest to us is redundant
        backward_index = valve_index - 2
        while backward_index >= 0 and 0 <= x1 - valve_locations[backward_index] <= MAX_DISTANCE_BETWEEN_VALVES:
          print(f"Removing {VALVE} due to backward redundancy where valve at {valve_locations[backward_index + 1]} is between two valves {valve_locations[backward_index]} and {x1} that are less than or equal to {MAX_DISTANCE_BETWEEN_VALVES} meters apart")
          del valve_locations[backward_index + 1]
          valve_index -= 1 # move index back one step due to removal
          backward_index = valve_index - 2

        # note that when dealing with forward redundant valves,
        # we must choose the valve within the range as the valve
        # we decide on removing or not
        forward_index = valve_index + 1
        while forward_index < len(valve_locations) and 0 <= valve_locations[forward_index] - x1 <= MAX_DISTANCE_BETWEEN_VALVES:
          print(f"Removing {VALVE} due to forward redundancy where valve at {valve_locations[forward_index - 1]} is between two valves {x1} and {valve_locations[forward_index]} that are less than or equal to {MAX_DISTANCE_BETWEEN_VALVES} meters apart")
          del valve_locations[forward_index - 1]


      # check if valve is within range
      while valve_index < len(valve_locations) and x1 <= valve_locations[valve_index] <= x2:
        valve_index += 1
      
      # finally, add the new pressure and cost for this section
      pressure += pressure_delta
      total_cost += section_cost
    
    # checks if we are at the end of pipeline (not terminated)
    if valve_index == len(valve_locations):
      # once the pipeline is built, we must find the pressure
      # at the final point in E5. This pressure will be the
      # last recorded pressure minus minor losses at outlet
      pressure -= get_minor_loss(KL_PIPE_OUTLET, P, V_AVG)
      print("PRESSURE AT OUTLET:", get_kilo_pascals(pressure), "kPa")
      pressure_data.append((x1, pressure))

      # get the total cost of the pipeline
      print(f"TOTAL COST: ${total_cost}")
    
    # update the minimum cost if this pipe type is
    if total_cost < minimum_cost:
      minimum_cost = total_cost
      minimum_cost_pipe_type = pipe_type
      minimum_cost_pressure_data = pressure_data
      minimum_cost_pump_station_locations = pump_station_locations
  
  print() # spacer for legibility
  print("===== RESULTS =====")
  print(f"Minimum Cost: ${minimum_cost}")
  print(f"Minimum Cost Pipe Type: {minimum_cost_pipe_type}")

  # add headers to export data
  PUMP_STATIONS_HEADER = ('x', 'elevation')
  PRESSURE_DATA_HEADER = ('x', 'pressure')
  minimum_cost_pump_station_locations.insert(0, PUMP_STATIONS_HEADER)
  minimum_cost_pressure_data.insert(0, PRESSURE_DATA_HEADER)
  
  # export data
  print() # spacer for legibility
  print("===== EXPORT DATA =====")
  export_data('Pump Stations', 'output_data/pump_stations.csv', minimum_cost_pump_station_locations)
  export_data('Pressure Data', 'output_data/pressure_data.csv', minimum_cost_pressure_data)   

# ========== SCRIPT END ==========








# ========== AUXILIARY START ==========

# ---------- ALL POSSIBLE REYNOLDS NUMBERS BY VARYING PIPE SIZE ----------
# this shows that no matter the pipe size, the flow will
# always remain turbulent as the reynold's number will
# ALWAYS be greater than 4000

# for pipe_type in PIPE_TYPES.keys():
#   ID = get_inner_diameter(pipe_type)
#   A = get_cross_sectional_area(ID)
#   V_AVG = get_average_velocity(V_RATE, A)
#   RE = get_reynolds_number(P, V_AVG, ID, MU)
#   print(RE)
  
# ========== AUXILIARY END ==========