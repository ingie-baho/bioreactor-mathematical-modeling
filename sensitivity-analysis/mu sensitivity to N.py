import numpy as np
import matplotlib.pyplot as plt

# Constants
mu_max = 2.8  # day^-1
Ks_N = 1.4  # g/m^3
P = 0.00693  # g/m^3
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

# Varying N from 0.247 to 2
N_range = np.linspace(0.247, 30, 100)
mu_values = [mu(N, P, Ci, I) for N in N_range]

# Plotting
plt.figure(figsize=(10, 5))
plt.plot(N_range, mu_values, 'b-')
plt.title('Growth Rate µ as a Function of Nutrient Concentration N', fontsize=18)  # Increase title font size if needed
plt.xlabel('Nutrient Concentration N (g/m^3)', fontsize=16)  # Set the fontsize for x-axis label
plt.ylabel('Growth Rate µ (day^-1)', fontsize=16)  # Set the fontsize for y-axis label
plt.grid(True)
plt.legend()
plt.show()