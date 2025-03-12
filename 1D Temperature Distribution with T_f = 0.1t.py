import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve


DT = 1.0  
Da = 1.0  
c = 1.0  
k = 1.0   
beta = 1.0  

L = 1 
Nx = 100  
Nt = 10000 
Dx = 2 * L / (Nx - 1)  
Dt = 3 / Nt 

Q = 0.1

x = np.linspace(-L, L, Nx)
t = np.linspace(0, 3, Nt)  

T = np.zeros((Nx, Nt))
A = np.ones((Nx, Nt))  

diagonal_T = (1 + 2 * DT * Dt / Dx**2) * np.ones(Nx)
off_diagonal_T = (-DT * Dt / Dx**2) * np.ones(Nx - 1)
D2_implicit_T = diags([diagonal_T, off_diagonal_T, off_diagonal_T], [0, -1, 1], shape=(Nx, Nx))

diagonal_A = (1 + 2 * Da * Dt / Dx**2) * np.ones(Nx)
off_diagonal_A = (-Da * Dt / Dx**2) * np.ones(Nx - 1)
D2_implicit_A = diags([diagonal_A, off_diagonal_A, off_diagonal_A], [0, -1, 1], shape=(Nx, Nx))

D2_implicit_T = D2_implicit_T.tolil()
D2_implicit_T[0, :] = 0
D2_implicit_T[-1, :] = 0
D2_implicit_T[0, 0] = 1
D2_implicit_T[-1, -1] = 1
D2_implicit_T = D2_implicit_T.tocsc()  

D2_implicit_A = D2_implicit_A.tolil()
D2_implicit_A[0, 0] = -1
D2_implicit_A[0, 1] = 1
D2_implicit_A[-1, -1] = -1
D2_implicit_A[-1, -2] = 1
D2_implicit_A = D2_implicit_A.tocsc()  

for n in range(Nt - 1):

    reaction_T = c * A[:, n] * np.exp(np.clip(beta * T[:, n], -50, 50)) * Dt  
   
    rhs_T = T[:, n] + reaction_T 
    T[:, n+1] = spsolve(D2_implicit_T, rhs_T)  

    T_f = Q*t[n+1]
    T[0, n+1] = T_f
    T[-1, n+1] = T_f

    reaction_A = -k * A[:, n] * np.exp(beta * T[:, n]) * Dt

    rhs_A = A[:, n] + reaction_A  
    A[:, n+1] = spsolve(D2_implicit_A, rhs_A) 

    A[0, n+1] = A[1, n+1] 
    A[-1, n+1] = A[-2, n+1] 

plt.figure(figsize=(8, 6))
plt.imshow(T.T, aspect='auto', extent=[x[0], x[-1], 0, 3], origin='lower', cmap='YlOrRd')
plt.colorbar(label='Temperature')
plt.xlabel(r"$x$", fontsize=14)
plt.ylabel(r"$t$", fontsize=14)
plt.show()

