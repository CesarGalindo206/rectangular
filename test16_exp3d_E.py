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

# --- 2. DEFINICIÓN DE LA MALLA (GRID) ---
puntos = 200
E_vals = np.linspace(0.001, 0.6, puntos)  
V0_vals = np.linspace(0.001, 0.6, puntos)  

# Creamos las matrices de coordenadas
E_grid, V0_grid = np.meshgrid(E_vals, V0_vals)

# Inicializamos la matriz de Transmitancia
T_grid = np.zeros((puntos, puntos))

# --- 3. CÁLCULO ITERATIVO (Mantenido igual) ---

for i in range(len(V0_vals)):
    V0 = V0_vals[i]
    for j in range(len(E_vals)):
        E = E_vals[j]
        
        K = np.sqrt(C1 * E)
        
        # --- LÓGICA CONDICIONAL ---
        if E < V0: # Efecto Túnel
            Q = np.sqrt(C2 * (V0 - E))
            if Q == 0: Q = 1e-15
            
            num = ((K * gamma)**2 + Q**2)**2
            den = 4 * (K**2) * (Q**2) * (gamma**2)
            term = (num / den) * np.sinh(Q * L)**2
            T_grid[i, j] = 1 / (1 + term)
            
        elif E > V0: # Propagación
            Q = np.sqrt(C2 * (E - V0))
            if Q == 0: Q = 1e-15
            
            num = ((K * gamma)**2 - Q**2)**2
            den = 4 * (K**2) * (Q**2) * (gamma**2)
            term = (num / den) * np.sin(Q * L)**2
            T_grid[i, j] = 1 / (1 + term)
            
        else: # Borde
            T_grid[i, j] = 0.5

# --- 4. GRAFICACIÓN MODIFICADA (E en eje Z) ---

fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

# CAMBIO DE EJES:
# X = V0_grid (Altura Barrera)
# Y = T_grid  (Transmitancia)
# Z = E_grid  (Energía - ahora es la altura)
surf = ax.plot_surface(V0_grid, T_grid, E_grid, cmap='twilight', 
                        linewidth=0, antialiased=False, alpha=0.8)

# Etiquetas
ax.set_xlabel('Altura Barrera $V_0$ (eV)', fontsize=11)
ax.set_ylabel('Transmitancia $T$', fontsize=11)
ax.set_zlabel('Energía $E$ (eV)', fontsize=11) 

ax.set_title(rf'Transmisión 3-Dimensional $V_0$,T,E' + '\n' + rf'$L={L}$ nm, $\gamma={gamma:.2f}$', fontsize=14)

# Ajustes de límites
ax.set_xlim(0, 0.6)
ax.set_ylim(0, 1.01)
ax.set_zlim(0, 0.6)

# Barra de color (ahora representa Energía E)
fig.colorbar(surf, shrink=0.5, aspect=10, label='Energía $E$ (eV)')

# --- LÍNEA DE REFERENCIA E=V0 ---
# Coordenadas:
# X (V0) = varía linealmente
# Y (T)  = constante en 0.5 (punto de cruce teórico)
# Z (E)  = varía linealmente (igual a V0)

line_vals = np.linspace(0.001, 0.6, 100)
T_line_vals = np.full_like(line_vals, 0.5) # Array lleno de 0.5

ax.plot(line_vals, T_line_vals, line_vals, 
        label='$E=V_0$ ($T=0.5$)', color='gray', linestyle='--', linewidth=3)

# Vista sugerida
ax.view_init(elev=30, azim=-110)

plt.legend()
plt.tight_layout()
plt.show()