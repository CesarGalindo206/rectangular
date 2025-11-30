import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

# --- 1. PARÁMETROS FÍSICOS ---
F_unidades = 26.2468  # [eV^-1 nm^-2]
m1_eff = 0.067  
m2_eff = 0.092  
L = 6.0         # Es (b-a) nm

C1 = m1_eff * F_unidades
C2 = m2_eff * F_unidades
gamma = m2_eff / m1_eff


puntos = 200
E_vals = np.linspace(0.001, 0.6, puntos)  
V0_vals = np.linspace(0.001, 0.6, puntos)  


E_grid, V0_grid = np.meshgrid(E_vals, V0_vals)


T_grid = np.zeros((puntos, puntos))


for i in range(len(V0_vals)):
    V0 = V0_vals[i]

    for j in range(len(E_vals)):
        E = E_vals[j]
        

        K = np.sqrt(C1 * E)
        

        if E < V0:
            # CASO 1: Efecto Túnel
            # Q es real
            Q = np.sqrt(C2 * (V0 - E))
            

            if Q == 0: Q = 1e-15
            
   
            num = ((K * gamma)**2 + Q**2)**2
            den = 4 * (K**2) * (Q**2) * (gamma**2)
            term = (num / den) * np.sinh(Q * L)**2
            
            T_grid[i, j] = 1 / (1 + term)
            
        elif E > V0:
            # CASO 2: Propagación 

            Q = np.sqrt(C2 * (E - V0))
            
            if Q == 0: Q = 1e-15
            
            num = ((K * gamma)**2 - Q**2)**2
            den = 4 * (K**2) * (Q**2) * (gamma**2)
            term = (num / den) * np.sin(Q * L)**2
            
            T_grid[i, j] = 1 / (1 + term)
            
        else:
            # CASO 3: E == V0 
            T_grid[i, j] = 0.5



fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')


surf = ax.plot_surface(E_grid, V0_grid, T_grid, cmap='inferno', 
                        linewidth=0, antialiased=False, alpha=0.9)

ax.set_xlabel('Energía $E$ (eV)', fontsize=11)
ax.set_ylabel('Altura Barrera $V_0$ (eV)', fontsize=11)
ax.set_zlabel('Transmitancia $T$', fontsize=11)
ax.set_title(f'Transmisión 3-Dimensional E,$V_0$,T \n$L={L}$ nm, $\\gamma={gamma:.2f}$', fontsize=14)


ax.set_zlim(0, 1.01)


fig.colorbar(surf, shrink=0.5, aspect=10, label='Probabilidad de Transmisión')


line_vals = np.linspace(0.001, 0.6, 100)
ax.plot(line_vals, line_vals, zs=0, zdir='z', label='$E=V_0$', color='gray', linestyle='--', linewidth=2)

ax.view_init(elev=30, azim=-165)

plt.tight_layout()
plt.show()