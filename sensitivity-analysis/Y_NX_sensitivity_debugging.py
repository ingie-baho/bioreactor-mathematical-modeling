import numpy as np
import matplotlib.pyplot as plt

# Constants
mu_max = 2.8  # day^-1
N0 = 5  # initial concentration of N, g/m^3
Ks_N = 1.4  # g/m^3
P0 = 10  # initial concentration of P, g/m^3
Ks_P = 0.06  # g/m^3
Ci = 6  # mol/m^3
Ks_Ci = 0.0295  # mol/m^3
K_I_Ci = 1000  # mol/m^3
I = 14  # umol/m^2 s
Ks_I = 50.6  # umol/m^2 s
K_I_I = 800  # umol/m^2 s
kd = 0.1  # day^-1
X0 = 1  # g/m^3

# Time setup
dt = 1 / 24  # time step in days (1 hour)
t_final = 14  # total time in days
time = np.arange(0, t_final, dt)

def mu(N, P, Ci, I):
    term_N = N / (N + Ks_N)
    term_P = P / (P + Ks_P)
    term_Ci = Ci / (Ci + Ks_Ci + (Ci ** 2 / K_I_Ci))
    term_I = I / (I + Ks_I + (I ** 2 / K_I_I))
    return mu_max * term_N * term_P * term_Ci * term_I

def rk4(X, N, P, t, dt, Y_NX):
    current_mu = mu(N, P, Ci, I)
    dX = current_mu * X - kd * X
    dN = -Y_NX * current_mu * X
    k1_X = dt * dX
    k1_N = dt * dN

    k2_X = dt * (current_mu * (X + 0.5 * k1_X) - kd * (X + 0.5 * k1_X))
    k2_N = dt * (-Y_NX * current_mu * (X + 0.5 * k1_X))

    k3_X = dt * (current_mu * (X + 0.5 * k2_X) - kd * (X + 0.5 * k2_X))
    k3_N = dt * (-Y_NX * current_mu * (X + 0.5 * k2_X))

    k4_X = dt * (current_mu * (X + k3_X) - kd * (X + k3_X))
    k4_N = dt * (-Y_NX * current_mu * (X + k3_X))

    X += (k1_X + 2 * k2_X + 2 * k3_X + k4_X) / 6
    N += (k1_N + 2 * k2_N + 2 * k3_N + k4_N) / 6
    print(f"μ: {current_mu:.4f}, N: {N:.4f}, dN: {dN:.4f}")
    return X, N

Y_NX_values = [10, 20]  # Specific Y_NX values for detailed output

for Y_NX in Y_NX_values:
    X = np.zeros(len(time))
    N = np.full(len(time), N0)
    X[0] = X0
    print(f"\nSimulation for Y_NX = {Y_NX}")
    for i in range(1, len(time)):
        # Update X and N using the rk4 method
        X[i], N[i] = rk4(X[i-1], N[i-1], P0, time[i-1], dt, Y_NX)
        # Print the updated values of N after the rk4 call

    plt.plot(time, N, label=f"Y_NX = {Y_NX}")

plt.title("Nitrogen Concentration Over Time")
plt.xlabel("Time (days)")
plt.ylabel("Nitrogen Concentration (g/m^3)")
plt.legend()
plt.show()
