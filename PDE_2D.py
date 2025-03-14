#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 19:37:27 2025

"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

D_T = 1.0  
c = 2.0   
a_0 = 1.0  
beta = 1.0  
T_f = 0.0  


nx, ny = 100, 100  
dx, dy = 2 / (nx - 1), 2 / (ny - 1) 
dt = 0.0001  
nt = 12000   

T = np.zeros((nx, ny)) 
T[:, 0] = T_f   
T[:, -1] = T_f  
T[0, :] = T_f   
T[-1, :] = T_f  

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

x = np.linspace(-1, 1, nx)
y = np.linspace(-1, 1, ny)
X, Y = np.meshgrid(x, y)

z_min = 0
z_max = 1.2
ax.set_zlim(z_min, z_max)

for n in range(nt):
    Tn = T.copy()  

    T[1:-1, 1:-1] = Tn[1:-1, 1:-1] + dt * (
        D_T * ((Tn[2:, 1:-1] - 2*Tn[1:-1, 1:-1] + Tn[:-2, 1:-1]) / dx**2 
               +(Tn[1:-1, 2:] - 2*Tn[1:-1, 1:-1] + Tn[1:-1, :-2]) / dy**2) 
        + c * a_0 * np.exp(beta * Tn[1:-1, 1:-1]))

    T[:, 0] = T_f   
    T[:, -1] = T_f  
    T[0, :] = T_f  
    T[-1, :] = T_f 
    
    T_inverted = T

    if n % 2000 == 0:
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(X, Y, T_inverted, cmap='RdYlBu_r', rstride=1, cstride=1, alpha=0.8)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlim(z_min, z_max)
        ax.set_zlabel('Temperature')
        ax.set_title(f'Time-Dependent Temperature Distribution (t = {n * dt:.4f})')
        plt.show()
        
        