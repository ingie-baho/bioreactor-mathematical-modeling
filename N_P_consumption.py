import numpy as np
import matplotlib.pyplot as plt

# Constants
mu_max = 2.8  # day^-1
N0 = 20  # initial concentration of N, g/m^3
Ks_N = 1.4  # g/m^3
P0 = 25  # initial concentration of P, g/m^3
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

# Time setup
dt = 1 / 24  # time step in days (1 hour)
t_final = 7  # total time in days
time = np.arange(0, t_final, dt)

# Growth rate function
def mu(N, P, Ci, I):
    term_N = N / (N + Ks_N)
    term_P = P / (P + Ks_P)
    term_Ci = Ci / (Ci + Ks_Ci + (Ci ** 2 / K_I_Ci))
    term_I = I / (I + Ks_I + (I ** 2 / K_I_I))
    return mu_max * term_N * term_P * term_Ci * term_I

# RK4 Integration
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

# Simulation
X = np.zeros(len(time))
N = np.full(len(time), N0)
P = np.full(len(time), P0)
mu_values = np.zeros(len(time))
X[0], N[0], P[0] = X0, N0, P0
mu_values[0] = mu(N[0], P[0], Ci, I)

for i in range(1, len(time)):
    X[i], N[i], P[i], mu_val = rk4(X[i - 1], N[i - 1], P[i - 1], time[i - 1], dt)
    mu_values[i] = mu_val
# Plotting using subplots
plt.figure(figsize=(12, 10))

# Plot for biomass concentration
plt.subplot(2, 1, 1)  # 2 rows, 1 column, 1st subplot
plt.plot(time, X / 1000, label='Biomass (X)', color='blue', marker = 'o')  # Convert g/m^3 to kg/m^3
plt.title('Biomass Concentration Over Time', fontsize=18)
plt.xlabel('Time (days)', fontsize=16)
plt.ylabel('Biomass (kg/m^3)', fontsize=16)  # Updated unit
plt.grid(True)
plt.xlim(left=0)  # Setting the x-axis to start from 0
plt.legend()

# Plot for nitrogen and phosphorus concentrations
plt.subplot(2, 1, 2)  # 2 rows, 1 column, 2nd subplot
plt.plot(time, N / 1000, label='Nitrogen (N)', color='green')  # Convert g/m^3 to kg/m^3
plt.plot(time, P / 1000, label='Phosphorus (P)', color='red')  # Convert g/m^3 to kg/m^3
plt.title('Nitrogen and Phosphorus Concentrations Over Time', fontsize=18)
plt.xlabel('Time (days)', fontsize=16)
plt.ylabel('Concentration (kg/m^3)', fontsize=16)  # Updated unit
plt.grid(True)
plt.xlim(left=0)  # Setting the x-axis to start from 0
plt.legend()

plt.tight_layout()
plt.show()
