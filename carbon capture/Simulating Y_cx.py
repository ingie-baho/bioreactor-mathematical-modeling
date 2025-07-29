import numpy as np
import matplotlib.pyplot as plt

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
M_net = -0.01

# Constants for X calculation
mu_max = 2.8  # day^-1
N0 = 20       # g/m^3
P0 = 25       # g/m^3
Ks_N = 1.4    # g/m^3
Ks_P = 0.06   # g/m^3
Ci = 6        # mol/m^3
Ks_Ci = 0.0295  # mol/m^3
K_I_Ci = 1000   # mol/m^3
I = 14          # umol/m^2 s
Ks_I = 50.6     # umol/m^2 s
K_I_I = 800     # umol/m^2 s
kd = 0.1        # day^-1
X0 = 1          # g/m^3
V = 0.005       # volume in m^3

# Time setup
dt = 1 / 24  # days
t_final = 7  # days
time = np.arange(0, t_final, dt)

# Arrays for simulation
X = np.zeros(len(time))
N = np.full(len(time), N0)
P = np.full(len(time), P0)

X[0] = X0

# Calculating M_c
# Calculations for M_a and M_r (mol/s)
M_a = K_La_air * (C_L_air - C1) * Vt
M_r = G * (y_in - (R*T*H / p_r) * (C1 + (p * y_in / R*T*H - C1) * np.exp(-K_La_CO2 * (1 - epsilon) * A * p * h / G*R*T*H)))
M_c = M_a + M_r - M_net

# Function to compute mu
def mu(N, P, Ci, I):
    term_N = N / (N + Ks_N)
    term_P = P / (P + Ks_P)
    term_Ci = Ci / (Ci + Ks_Ci + (Ci ** 2 / K_I_Ci))
    term_I = I / (I + Ks_I + (I ** 2 / K_I_I))
    return mu_max * term_N * term_P * term_Ci * term_I

# Run simulation using RK4
for i in range(1, len(time)):
    # RK4 step to find X, N, P
    current_mu = mu(N[i-1], P[i-1], Ci, I)
    dX = current_mu * X[i-1] - kd * X[i-1]
    dN = -0.1 * current_mu * X[i-1]
    dP = -0.03 * current_mu * X[i-1]
    X[i] = X[i-1] + dt * dX
    N[i] = N[i-1] + dt * dN
    P[i] = P[i-1] + dt * dP

# Calculate Y using the formula given
dX_dt = np.gradient(X, dt)  # numerical derivative of X with respect to time
Y = (1 / M_c) * dX_dt * V

# Plotting results
plt.figure(figsize=(12, 8))
plt.subplot(2, 1, 1)
plt.plot(time, X, label='Biomass (X)', color='blue')
plt.title('Biomass Concentration Over Time ', fontsize=18)
plt.xlabel('Time (days)', fontsize=16)
plt.ylabel('Biomass (g/m³)', fontsize=16)
plt.grid(True)
plt.legend()

plt.subplot(2, 1, 2)
plt.plot(time, Y, label='Yield (Y)', color='green')
plt.title('Yield Over Time', fontsize=18)
plt.xlabel('Time (days)', fontsize=16)
plt.ylabel('Yield (Y) [g carbon/g biomass]', fontsize=16)
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()
