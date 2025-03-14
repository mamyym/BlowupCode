#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 10 23:42:51 2025

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jn
from mpl_toolkits.mplot3d import Axes3D

lambda_val = 2
C = 1
Tf = 1

theta = np.linspace(0, 2*np.pi, 100)  
r = np.linspace(0, 1, 50)

R, Theta = np.meshgrid(r, theta)
X = R * np.cos(Theta)  
Y = R * np.sin(Theta)

A = (Tf - (C / lambda_val)) / jn(0, np.sqrt(lambda_val))

T = A * jn(0, np.sqrt(lambda_val) * R) + (C / lambda_val)

fig = plt.figure(figsize=(7, 7))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, T, cmap="YlOrRd", alpha=0.8)

ax.set_xlabel("x", fontsize=14)
ax.set_ylabel("y", fontsize=14)
ax.set_zlabel("Temperature", labelpad=15, fontsize=14)
plt.show()  