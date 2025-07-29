import numpy as np
import matplotlib.pyplot as plt

# Constants
mu_max = 2.8  # day^-1
N = 30  # g/m^3
Ks_N = 1.4  # g/m^3
P = 2  # g/m^3
Ks_P = 0.06  # g/m^3
Ci = 6  # mol/m^3
Ks_Ci = 0.0295  # mol/m^3
K_I_Ci = 1000  # mol/m^3
I = 200  # umol/m^2 s
Ks_I = 50.6  # umol/m^2 s
K_I_I = 800  # umol/m^2 s
kd = 0.01  # day^-1
X0 = 1  # g/m^3

# Time setup
dt = 1/24  # time step in days (1 hour)
t_final = 14  # total time in days
time = np.arange(0, t_final, dt)

# Growth rate function mu
def mu(N, P, Ci, I):
    term_N = N / (N + Ks_N)
    term_P = P / (P + Ks_P)
    term_Ci = Ci / (Ci + Ks_Ci + (Ci**2 / K_I_Ci))
    term_I = I / (I + Ks_I + (I**2 / K_I_I))
    return mu_max * term_N * term_P * term_Ci * term_I

# RK4 Integration
def rk4(X, t, dt):
    current_mu = mu(N, P, Ci, I)
    k1 = dt * (current_mu * X - kd * X)
    k2 = dt * (current_mu * (X + 0.5 * k1) - kd * (X + 0.5 * k1))
    k3 = dt * (current_mu * (X + 0.5 * k2) - kd * (X + 0.5 * k2))
    k4 = dt * (current_mu * (X + k3) - kd * (X + k3))
    return X + (k1 + 2 * k2 + 2 * k3 + k4) / 6, current_mu

# Simulation
X = np.zeros(len(time))
mu_values = np.zeros(len(time))
X[0] = X0
mu_values[0] = mu(N, P, Ci, I)

for i in range(1, len(time)):
    X[i], mu_val = rk4(X[i-1], time[i-1], dt)
    mu_values[i] = mu_val

# Final mu value
final_mu = mu_values[-1]
comparison_result = "greater" if final_mu > kd else "less"

# Plotting
plt.figure(figsize=(10, 5))
plt.plot(time, X, label='Biomass X(t)')
plt.title('Biomass Concentration Over Time')
plt.xlabel('Time (days)')
plt.ylabel('Biomass (g/m^3)')
plt.grid(True)
plt.annotate(f'Final µ = {final_mu:.4f} (µ is {comparison_result} than k_d = {kd})',
             xy=(0.5, 0.1), xycoords='axes fraction',
             fontsize=12, bbox=dict(boxstyle="round,pad=0.3", edgecolor='orange', facecolor='white'))
plt.show()
