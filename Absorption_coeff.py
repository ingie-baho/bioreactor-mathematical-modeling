import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from math import log


# Constants
mu_max = 2.8  # day^-1
N0 = 5  # initial concentration of N, g/m^3
P0 = 5  # initial concentration of P, g/m^3
Ks_N = 1.4  # g/m^3
Ks_P = 0.06  # g/m^3
Ci = 6  # mol/m^3
Ks_Ci = 0.0295  # mol/m^3
K_I_Ci = 1000  # mol/m^3
I = 14  # umol/m^2 s
Ks_I = 50.6  # umol/m^2 s
K_I_I = 800  # umol/m^2 s
kd = 0.1  # day^-1
X0 = 1  # g/m^3
Y_NX = 0.1  # Yield of N on X
Y_PX = 0.03  # Yield of P on X
d = 0.004  # path length in meters
a = 30

# Time setup
dt = 1 / 24  # time step in days (1 hour)
t_final = 2  # total time in days
time = np.arange(0, t_final, dt)

# Import OD data
od_data = pd.read_csv('OD_log.txt', header=None, names=['timestamp', 'OD'], parse_dates=[0], infer_datetime_format=True)

# Ensure correct formatting if automatic parsing fails
od_data['timestamp'] = pd.to_datetime(od_data['timestamp'], format='%Y-%m-%d %H:%M:%S.%f')

# Calculate time in days relative to the start
od_data['time_days'] = (od_data['timestamp'] - od_data['timestamp'].iloc[0]).dt.total_seconds() / 86400

od_data['X'] = od_data['OD'] * log(10) / (d * a)

# Function definitions
def mu(N, P, Ci, I):
    term_N = N / (N + Ks_N)
    term_P = P / (P + Ks_P)
    term_Ci = Ci / (Ci + Ks_Ci + (Ci ** 2 / K_I_Ci))
    term_I = I / (I + Ks_I + (I ** 2 / K_I_I))
    return mu_max * term_N * term_P * term_Ci * term_I

def rk4(X, N, P, t, dt):
    current_mu = mu(N, P, Ci, I)
    dX = current_mu * X - kd * X
    dN = -Y_NX * current_mu * X
    dP = -Y_PX * current_mu * X
    k1_X = dt * dX
    k1_N = dt * dN
    k1_P = dt * dP

    k2_X = dt * (current_mu * (X + 0.5 * k1_X) - kd * (X + 0.5 * k1_X))
    k2_N = dt * (-Y_NX * current_mu * (X + 0.5 * k1_X))
    k2_P = dt * (-Y_PX * current_mu * (X + 0.5 * k1_X))

    k3_X = dt * (current_mu * (X + 0.5 * k2_X) - kd * (X + 0.5 * k2_X))
    k3_N = dt * (-Y_NX * current_mu * (X + 0.5 * k2_X))
    k3_P = dt * (-Y_PX * current_mu * (X + 0.5 * k2_X))

    k4_X = dt * (current_mu * (X + k3_X) - kd * (X + k3_X))
    k4_N = dt * (-Y_NX * current_mu * (X + k3_X))
    k4_P = dt * (-Y_PX * current_mu * (X + k3_X))

    X += (k1_X + 2 * k2_X + 2 * k3_X + k4_X) / 6
    N += (k1_N + 2 * k2_N + 2 * k3_N + k4_N) / 6
    P += (k1_P + 2 * k2_P + 2 * k3_P + k4_P) / 6
    return X, N, P, current_mu

# Simulation loop
X = np.zeros(len(time))
N = np.full(len(time), N0)
P = np.full(len(time), P0)
X[0], N[0], P[0] = X0, N0, P0

for i in range(1, len(time)):
    X[i], N[i], P[i], _ = rk4(X[i - 1], N[i - 1], P[i - 1], time[i - 1], dt)

# Interpolate biomass concentration at OD measurement times
od_data['X_interp'] = np.interp(od_data['time_days'], time, X)
od_data['absorbance'] = od_data['OD'] * np.log(10) / (od_data['X_interp'] * d)

# Plotting absorbance over time
plt.figure(figsize=(12, 6))
plt.plot(od_data['time_days'], od_data['absorbance'], label='Absorbance', color='purple')
plt.title('Absorbance Over Time', fontsize=18)
plt.xlabel('Time (days)', fontsize=16)
plt.ylabel('Absorbance [m^2/g]', fontsize=16)
plt.grid(True)
plt.legend()
plt.show()

# Plotting
plt.figure(figsize=(10, 5))
plt.plot(od_data['time_days'], od_data['X'], label='Biomass (X) calculated from OD', color='blue')
plt.title('Biomass Concentration Over Time Calculated from OD', fontsize=18)
plt.xlabel('Time (days)', fontsize=16)
plt.ylabel('Biomass Concentration X (g/m^3)', fontsize=16)
plt.grid(True)
plt.legend()
plt.show()

# Ensure the rest of your plotting code is intact as needed
