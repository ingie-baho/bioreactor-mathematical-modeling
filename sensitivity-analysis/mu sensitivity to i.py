import numpy as np
import matplotlib.pyplot as plt

# Constants
mu_max = 2.8  # day^-1
N = 0.247  # g/m^3
Ks_N = 1.4  # g/m^3
P = 0.00693  # g/m^3
Ks_P = 0.06  # g/m^3
Ci = 2  # mol/m^3
Ks_Ci = 0.0295  # mol/m^3
K_I_Ci = 1000  # mol/m^3
Ks_I = 50.6  # umol/m^2 s
K_I_I = 800  # umol/m^2 s

# Growth rate function mu
def mu(N, P, Ci, I):
    term_N = N / (N + Ks_N)
    term_P = P / (P + Ks_P)
    term_Ci = Ci / (Ci + Ks_Ci + (Ci**2 / K_I_Ci))
    term_I = I / (I + Ks_I + (I**2 / K_I_I))
    return mu_max * term_N * term_P * term_Ci * term_I

# Varying I from 14 to 500
I_range = np.linspace(14, 500, 100)
mu_values = [mu(N, P, Ci, I) for I in I_range]

# Plotting
plt.figure(figsize=(10, 5))
plt.plot(I_range, mu_values, 'm-', label='Growth Rate µ(I)')
plt.title('Growth Rate µ as a Function of Light Intensity I', fontsize=18)
plt.xlabel('Light Intensity I (μmol/m^2 s)', fontsize=16)
plt.ylabel('Growth Rate µ (day^-1)', fontsize=16)
plt.grid(True)
plt.legend()
plt.show()
