import numpy as np
import matplotlib.pyplot as plt

# Constants
mu_max = 2.8  # day^-1
N = 0.247  # g/m^3
Ks_N = 1.4  # g/m^3
Ks_P = 0.06  # g/m^3
Ci = 2  # mol/m^3
Ks_Ci = 0.0295  # mol/m^3
K_I_Ci = 1000  # mol/m^3
I = 14  # umol/m^2 s
Ks_I = 50.6  # umol/m^2 s
K_I_I = 800  # umol/m^2 s

# Growth rate function mu
def mu(N, P, Ci, I):
    term_N = N / (N + Ks_N)
    term_P = P / (P + Ks_P)
    term_Ci = Ci / (Ci + Ks_Ci + (Ci**2 / K_I_Ci))
    term_I = I / (I + Ks_I + (I**2 / K_I_I))
    return mu_max * term_N * term_P * term_Ci * term_I

# Varying P from 0.00693 to 5
P_range = np.linspace(0.00693, 5, 100)
mu_values = [mu(N, P, Ci, I) for P in P_range]

# Plotting
plt.figure(figsize=(10, 5))
plt.plot(P_range, mu_values, 'g-', label='Growth Rate µ(P)')
plt.title('Growth Rate µ as a Function of Phosphorus Concentration P')
plt.xlabel('Phosphorus Concentration P (g/m^3)')
plt.ylabel('Growth Rate µ (day^-1)')
plt.grid(True)
plt.legend()
plt.show()
