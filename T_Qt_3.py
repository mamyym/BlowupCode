#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 7 14:41:23 2025

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
Dt = 218.8 / Nt 

x = np.linspace(-L, L, Nx)
t = np.linspace(0, 218.8, Nt) 

Q_values = [0.1, 0.01]
Tf_constant = 3 

temperature_evolutions = {}

for Q in Q_values:
    T = np.zeros((Nx, Nt))

    diagonal = (1 + 2 * DT * Dt / Dx**2) * np.ones(Nx)
    off_diagonal = (-DT * Dt / Dx**2) * np.ones(Nx - 1)
    D2_implicit = diags([diagonal, off_diagonal, off_diagonal], [0, -1, 1], shape=(Nx, Nx))

    D2_implicit = D2_implicit.tolil()
    D2_implicit[0, :] = 0
    D2_implicit[-1, :] = 0
    D2_implicit[0, 0] = 1
    D2_implicit[-1, -1] = 1
    D2_implicit = D2_implicit.tocsc() 

    for n in range(Nt - 1):
       
        reaction = c * a * np.exp(beta * T[:, n]) * Dt

        rhs = T[:, n] + reaction 
        T[:, n+1] = spsolve(D2_implicit, rhs)  
        Tf = Q * t[n+1]
        T[0, n+1] = Tf
        T[-1, n+1] = Tf

    center_index = Nx // 2 
    temperature_evolutions[f"Q = {Q}"] = T[center_index, :]

T_const = np.zeros((Nx, Nt))
D2_implicit = diags([diagonal, off_diagonal, off_diagonal], [0, -1, 1], shape=(Nx, Nx))
D2_implicit = D2_implicit.tolil()
D2_implicit[0, :] = 0
D2_implicit[-1, :] = 0
D2_implicit[0, 0] = 1
D2_implicit[-1, -1] = 1
D2_implicit = D2_implicit.tocsc()

for n in range(Nt - 1):
    reaction = c * a * np.exp(beta * T_const[:, n]) * Dt
    rhs = T_const[:, n] + reaction
    T_const[:, n+1] = spsolve(D2_implicit, rhs)
    T_const[0, n+1] = Tf_constant
    T_const[-1, n+1] = Tf_constant


temperature_evolutions[f"T_f = {Tf_constant}"] = T_const[center_index, :]

plt.figure(figsize=(8, 5))
for label, T_center in temperature_evolutions.items():
    plt.plot(t, T_center, label=label)
plt.xlabel("$t$", fontsize=14)
plt.ylabel("$T(0,t)$", fontsize=14)
plt.xscale("log") 
plt.yscale("log") 
plt.ylim(1e-3, 1e5)  
plt.legend(fontsize=14)
plt.grid()
plt.show()
