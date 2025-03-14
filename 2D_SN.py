#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 6 00:58:48 2025

"""

import numpy as np
import matplotlib.pyplot as plt

D_T = 1.0 
c = 1.0    
a_0 = 1.0  
beta = 0.1 
T_f = 1.0  

nx, ny = 50, 50  
dx, dy = 2 / (nx - 1), 2 / (ny - 1)  
dt = 0.0001  
nt = 1000    

T = np.zeros((nx, ny))  
T[:, 0] = T_f   
T[:, -1] = T_f  
T[0, :] = T_f   
T[-1, :] = T_f  

for n in range(nt):
    Tn = T.copy() 

    T[1:-1, 1:-1] = Tn[1:-1, 1:-1] + dt * (
        D_T * (
            (Tn[2:, 1:-1] - 2*Tn[1:-1, 1:-1] + Tn[:-2, 1:-1]) / dx**2 + 
            (Tn[1:-1, 2:] - 2*Tn[1:-1, 1:-1] + Tn[1:-1, :-2]) / dy**2
        ) + c * a_0 * np.exp(beta * Tn[1:-1, 1:-1])
    )

    T[:, 0] = T_f   
    T[:, -1] = T_f  
    T[0, :] = T_f   
    T[-1, :] = T_f  

x = np.linspace(-1, 1, nx)
y = np.linspace(-1, 1, ny)
X, Y = np.meshgrid(x, y)

plt.figure(figsize=(8, 6))
plt.contourf(X, Y, T, levels=50, cmap='YlOrRd_r')
plt.colorbar(label='Temperature')
plt.xlabel('x')
plt.ylabel('y')
plt.show()