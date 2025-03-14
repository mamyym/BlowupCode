"""
Created on Thu Feb 20 19:37:27 2025

"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


D_T = 1.0   
c   = 1.0   
a_0 = 1.0  
beta= 1.0   
T_f = 3.0   
#3 - 23400(234800)

nx, ny = 100, 100
dx, dy = 2/(nx-1), 2/(ny-1)
dt     = 0.00001
nt     =23400

T = np.zeros((nx, ny)) 

T[:, 0]  = T_f
T[:, -1] = T_f
T[0, :]  = T_f
T[-1, :] = T_f

fractions = [0, 1/4, 1/2, 3/4, 99/100, 1.0]
snap_steps = [nt * frac for frac in fractions]
snapshots  = {}  

snapshots[0] = T.copy()

for n in range(1, nt+1):
    Tn = T.copy()

    T[1:-1, 1:-1] = Tn[1:-1, 1:-1] + dt * (
        D_T * (
            (Tn[2:,   1:-1] - 2*Tn[1:-1, 1:-1] + Tn[:-2, 1:-1]) / dx**2 +
            (Tn[1:-1, 2:]   - 2*Tn[1:-1, 1:-1] + Tn[1:-1, :-2]) / dy**2
        )
        + c * a_0 * np.exp(beta * Tn[1:-1, 1:-1])
    )
    

    T[:,  0] = T_f
    T[:, -1] = T_f
    T[0,  :] = T_f
    T[-1, :] = T_f
    

    if n in snap_steps:
        snapshots[n] = T.copy()

x = np.linspace(-1, 1, nx)
y = np.linspace(-1, 1, ny)
X, Y = np.meshgrid(x, y)

z_min = 0
z_max = 9

fig = plt.figure(figsize=(18, 10))

axes = []
for i in range(6):
    ax = fig.add_subplot(2, 3, i+1, projection='3d')
    axes.append(ax)

for i, n_step in enumerate(snap_steps):
    ax = axes[i]

    T_plot = snapshots[n_step]

    surf = ax.plot_surface(
        X, Y, T_plot, 
        cmap="YlOrRd", 
        rstride=1, 
        cstride=1, 
        alpha=0.9, 
        vmin=z_min, 
        vmax=z_max
    )
    
    ax.set_title(f"t = {n_step * dt:.3f}", fontsize=16)
    ax.set_xlabel("x", fontsize=14)
    ax.set_ylabel("y", fontsize=14)
    ax.set_zlabel("Temperature", fontsize=14, labelpad=20)
    ax.tick_params(axis='both', labelsize=14)  
    ax.tick_params(axis='z', labelsize=14, pad=10)    
    ax.set_zlim(z_min, z_max)

cax = fig.add_axes([1, 0.1, 0.01, 0.8])

cbar = fig.colorbar(
    surf, 
    cax=cax,
    shrink=0.7, 
    aspect=20
)
cbar.set_label("Temperature", fontsize=16)
cbar.ax.tick_params(labelsize=12)

plt.tight_layout()
plt.show()
