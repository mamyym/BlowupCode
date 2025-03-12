import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt


beta = 1
C = 0.5
D_T = 5

def Tf(T0):
    return T0 - (2 / beta) * np.log(np.cosh(np.sqrt(beta * C / (2 * D_T)) * np.exp(beta * T0 / 2)))

def dTf_dT0(T0):
    dT0 = 1e-5  
    return (Tf(T0 + dT0) - Tf(T0)) / dT0  

T0_crit = opt.fsolve(dTf_dT0, 3)[0]  
Tf_crit = Tf(T0_crit)  

T0_values = np.linspace(-2, 8, 400)
Tf_values = Tf(T0_values)

plt.figure(figsize=(8,6))

plt.plot(T0_values, Tf_values, color='black', label=r'$y = T_0 - \frac{2}{\beta} \ln \left[ \cosh \left( \sqrt{\frac{\beta c_{new}}{2D_T}} e^{\frac{\beta T_0}{2}} \right) \right]$')

slope = dTf_dT0(T0_crit)
T0_tangent = np.linspace(T0_crit - 2, T0_crit + 2, 100)
Tf_tangent = Tf_crit + slope * (T0_tangent - T0_crit)

plt.axhline(y=Tf_crit, color='red', linestyle='--', label=r'Critical Line: $y = T_c$')
plt.axhline(y=3, linestyle='--', color='green')
plt.axhline(y=1, linestyle='--', color='b')

plt.xlim(-2,6.5)
plt.ylim(-2,4.5)

plt.xlabel(r"$T_0$", fontsize=16)
plt.ylabel(r"$y$", fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=14)
plt.grid()


plt.show()