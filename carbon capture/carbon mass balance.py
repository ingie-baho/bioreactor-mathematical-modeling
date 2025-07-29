import numpy as np

# Constants and Variables for M_a

C_L_air = 0.012  # (mol/m³) eqlbm concentration of CO2 in water at 25°C
Vt = 0.005        # volume of the tank, in m^3
C1 = 0.15  #liquid phase concentration of CO2 (mol/m3)

K_La_air = 1.24e-4  # from paper, in 1/s (1.24e-4)

# Constants and Variables for M_r
G = 0.00123         # overall mass transfer coefficient [mol/sec]
R = 8.206e-5  # Idea gas constant [atm·m^3·mol^-1·K^-1]
T = 298.15      #room temp [kelvin]
p = 1         # total pressure surrounding reactor [atm]
y_in = 0.000420     #mole fraction of CO2 at the sparger [unitless] (0.000420)
h = 0.2032       #riser depth [m] (converted form 8 inches)
G = 0.00123         # gas flow rate at the sparger [mol/sec]
A = 0.018           # cross sectional area of the bioreactor tank [m^2]
epsilon = 0.006   # fraction of reactor's volume that's occupied by gas [unitless]
p_r = 1      # atmospheric of the sparged gas in the reactor [atm]

H = 0.034       #henry's constant for CO2 [unitless]
K_La_CO2 = 1.24e-4 # mass transfer coefficient for CO2, in 1/s
# Make up a value for M_net (mol/s)
M_net = -1.0


# Calculations for M_a and M_r (mol/s)
M_a = K_La_air * (C_L_air - C1) * Vt
M_r = G * (y_in - (R*T*H / p_r) * (C1 + (p * y_in / R*T*H - C1) * np.exp(-K_La_CO2 * (1 - epsilon) * A * p * h / G*R*T*H)))



# Calculate M_c
M_c = M_a + M_r - M_net

print(f"M_a: {M_a}")
print(f"M_r: {M_r}")
print(f"M_c: {M_c}")
