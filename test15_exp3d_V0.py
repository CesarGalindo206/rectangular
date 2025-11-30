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
        
        if E < V0: # Túnel
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

# --- 4. GRAFICACIÓN MODIFICADA (V0 en eje Z) ---

fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

# CAMBIO CLAVE AQUÍ:
# X = E_grid
# Y = T_grid (Ahora la Transmitancia está en el plano base/profundidad)
# Z = V0_grid (Ahora la Altura de Barrera es la altura del gráfico)
surf = ax.plot_surface(E_grid, T_grid, V0_grid, cmap='cividis', 
                        linewidth=0, antialiased=False, alpha=0.8)

# Etiquetas actualizadas
ax.set_xlabel('Energía $E$ (eV)', fontsize=11)
ax.set_ylabel('Transmitancia $T$', fontsize=11) 
ax.set_zlabel('Altura Barrera $V_0$ (eV)', fontsize=11) # Z es V0 ahora

ax.set_title(rf'Transmisión 3-Dimensional E,T,$V_0$' + '\n' + rf'$L={L}$ nm, $\gamma={gamma:.2f}$', fontsize=14)

# Ajustamos límites
ax.set_ylim(0, 1.01)       # T va de 0 a 1
ax.set_zlim(0, 0.6)        # V0 va de 0 a 0.6

# Barra de color (ahora representa V0, ya que es la altura Z)
fig.colorbar(surf, shrink=0.5, aspect=10, label='Altura de Barrera $V_0$ (eV)')

# --- DIBUJAR LA LÍNEA E=V0 EN ESTE NUEVO ESPACIO ---
# En la línea donde E = V0, sabemos que T = 0.5 siempre.
# Por lo tanto, las coordenadas de la línea roja son:
# X = line_vals (Energía)
# Y = 0.5 (Transmitancia constante en el borde)
# Z = line_vals (Porque V0 = E)

line_vals = np.linspace(0.001, 0.6, 100)
T_line_vals = np.full_like(line_vals, 0.5) # Array lleno de 0.5

ax.plot(line_vals, T_line_vals, line_vals, 
        label='$E=V_0$ (donde $T=0.5$)', color='gray', linestyle='--', linewidth=3)

# Vista sugerida para apreciar la estructura doblada
ax.view_init(elev=25, azim=-145)

plt.legend()
plt.tight_layout()
plt.show()