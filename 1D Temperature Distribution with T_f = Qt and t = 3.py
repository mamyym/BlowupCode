import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve


DT = 5.0  
c = 1.0   
a = 0.5 
beta = 1.0 

L = 1  
Nx = 100  
Nt = 10000  
Dx = 2 * L / (Nx - 1)  
Dt = 3 / Nt  

x = np.linspace(-L, L, Nx)
t = np.linspace(0, 3, Nt) 

Q_values = [0.1, 0.01, 0.001]

fig, axs = plt.subplots(1, 3, figsize=(15, 5), sharey=True)

im_list = []

for i, Q in enumerate(Q_values):
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
    
    vmin, vmax = np.min(T), np.max(T)

    im = axs[i].imshow(T.T, aspect='auto', extent=[x[0], x[-1], 0, 3], origin='lower', cmap='YlOrRd', vmin=vmin, vmax=vmax)
    im_list.append(im)
   
    axs[i].set_title(f"Q = {Q}", fontsize=14)

    axs[i].set_xlabel(r"$x$", fontsize=14)
    if i == 0:
        axs[i].set_ylabel(r"$t$", fontsize=14)  

plt.subplots_adjust(right=0.85)

cbar_ax = fig.add_axes([0.88, 0.15, 0.02, 0.7])  
cbar = fig.colorbar(im_list[-1], cax=cbar_ax)
cbar.set_label("Temperature", fontsize=14)

plt.savefig("temperature_evolution_Q_values_adjusted.png", dpi=300)


plt.show()
plt.show()