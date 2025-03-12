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

T0_vals = np.linspace(-2, 8, 400)
Tf_vals = np.linspace(-2, 3, 400)
T_0, T_f = np.meshgrid(T0_vals, Tf_vals)

def f(T_0, T_f):
    return T_f - T_0 + (2 / beta) * np.log(np.cosh(np.sqrt(beta * C / (2 * D_T)) * np.exp(beta * T_0 / 2)))

plt.figure(figsize=(8,6))

plt.contour(T_f, T_0, f(T_0, T_f), levels=[0], colors='black', linewidths=2)

plt.plot(Tf_crit, T0_crit, 'ro', markersize=8, label="Critical Point")

slope = dTf_dT0(T0_crit)
T0_tangent = np.linspace(T0_crit - 2, T0_crit + 2, 100)
Tf_tangent = Tf_crit + slope * (T0_tangent - T0_crit)

plt.axvline(x=Tf_crit, color='b', linestyle='--', label=r'Critical Line: $T_f = T_c$')

plt.xlabel(r"$T_f$", fontsize=16)
plt.ylabel(r"$T_0$", fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=14)
plt.grid()


plt.show()