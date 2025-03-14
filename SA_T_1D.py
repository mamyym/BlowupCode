#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  9 18:24:51 2025

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve


DT = 1.0     
c = 1.0    
a = 1.0      
beta = 1.0  

L = 1       
Nx = 100     
Nt = 100000  
Dx = 2 * L / (Nx - 1)
total_time = 0.328603
# DT: 0.5 - 0.351213; 1 - 0.328603; 2 - 0.2639221; 3 - 0.2217167
# c: 0.1 - 1.48066; 0.5 - 0.527843; 1 - 0.328603; 2 - 0.175606; 3 - 0.1174755
# beta: 0.5 - 1.0885; 2 - 0.021225
Dt = total_time / Nt

x = np.linspace(-L, L, Nx)
t = np.linspace(0, total_time, Nt)
T = np.zeros((Nx, Nt))

diag_main = (1 + 2 * DT * Dt / Dx**2) * np.ones(Nx)
diag_off  = (-DT * Dt / Dx**2) * np.ones(Nx-1)
D2_implicit = diags([diag_main, diag_off, diag_off], [0, -1, 1], shape=(Nx, Nx)).tolil()

D2_implicit[0, :] = 0
D2_implicit[-1, :] = 0
D2_implicit[0, 0] = 1
D2_implicit[-1, -1] = 1
D2_implicit = D2_implicit.tocsc()

for n in range(Nt - 1):

    reaction = c * a * np.exp(beta * T[:, n]) * Dt
    rhs = T[:, n] + reaction
    rhs[0] = 3
    rhs[-1] = 3
    T[:, n+1] = spsolve(D2_implicit, rhs)

fractions = [1/4, 1/2, 9/10, 1.0]
time_points = [t[-1] * frac for frac in fractions]


plt.figure(figsize=(8, 6))
for tp in time_points:
    idx = np.argmin(np.abs(t - tp))
    plt.plot(x, T[:, idx], label=f"t = {t[idx]:.3f}")
plt.xlabel("x", fontsize=16)
plt.ylabel("Temperature", fontsize=16)
plt.legend()
plt.grid(True)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.tight_layout()
plt.legend(fontsize=16)
plt.show()

# Figure 2: Temperature evolution at the center (x = 0) over time
center_index = Nx // 2
plt.figure(figsize=(8, 6))
plt.plot(t, T[center_index, :], color='blue')
plt.xlabel("t", fontsize=16)
plt.ylabel("Temperature", fontsize=16)
plt.grid(True)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.tight_layout()
plt.show()
