import numpy as np
import matplotlib.pyplot as plt


def f(T, lambd, beta=1):
    return lambd - np.exp(beta * T) / (beta * T)

T_vals = np.linspace(0, 3, 400)
lambda_vals = np.linspace(2, 5, 400)
T, lambd = np.meshgrid(T_vals, lambda_vals)

T_crit = 1  
lambda_crit = np.exp(1) / T_crit  

plt.figure(figsize=(8,6))
plt.contour(lambd, T, f(T, lambd), levels=[0], colors='black', linewidths=2, label="Bifurcation Curve")

plt.plot(lambda_crit, T_crit, 'ro', markersize=8, label="Critical Point")

plt.axvline(x=lambda_crit, color='b', linestyle='--', label=r'Critical Line: $\lambda = \lambda_c$')

plt.xlabel("λ", fontsize=16)  
plt.ylabel("T", fontsize=16)
plt.xticks(fontsize=14)  
plt.yticks(fontsize=14)
plt.legend(fontsize=14)  


plt.grid()

plt.show()