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
Nt = 1000  
Dx = 2 * L / (Nx - 1)  
Dt = 1 / Nt 

x = np.linspace(-L, L, Nx)
t = np.linspace(0, 1, Nt) 

temperature_evolutions = {}

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

    T[0, n+1] = -1
    T[-1, n+1] = -1

center_index = Nx // 2 
temperature_evolutions[0] = T[center_index, :]

fig, axs = plt.subplots(1, 2, figsize=(10, 4))

axs[1].plot(t, temperature_evolutions[0], color='blue', label=r"$T(0,t)$")
axs[1].set_title("(b) Temperature Evolution at $x=0$ with $T_f = -1$ and $t = 1$")
axs[1].set_xlabel("$t$", fontsize=14)
axs[1].set_ylabel("$T(0,t)$", fontsize=14)
axs[1].grid()


im = axs[0].imshow(T.T, aspect='auto', extent=[x[0], x[-1], 0, 1], origin='lower', cmap='YlOrRd')
axs[0].set_title("(a) 1D Temperature Distribution with $T_f = -1$ and $t = 1$")
axs[0].set_xlabel("$x$", fontsize=14)
axs[0].set_ylabel("$t$", fontsize=14)
fig.colorbar(im, ax=axs[0], label="Temperature")

plt.tight_layout()

plt.savefig("two_plots.png", dpi=300)


plt.show()

