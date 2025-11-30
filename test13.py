import numpy as np
import matplotlib.pyplot as plt

F_unidades = 26.2468  # Unidades en heteroestucturas [eV^-1 nm^-2] (2m0/hbar^2)

# Masas efectivas
m1_eff = 0.067  
m2_eff = 0.092  
gamma = m2_eff / m1_eff

# Factores C = 2m* / hbar^2
C1 = m1_eff * F_unidades
C2 = m2_eff * F_unidades

V0 = 0.23      # eV
L = 10       # Es (b-a) en nm 
Energias = np.linspace(0.001, 0.6, 1000)

T_lista = []
R_lista = []
Suma_lista = []

for E in Energias:
    
    K = np.sqrt(C1 * E)
    
    term_comun = 0.0
    
    if E < V0:
        
        Q = np.sqrt(C2 * (V0 - E))
        
        
        numerador = ((K * gamma)**2 + Q**2)**2
        denominador = 4 * (K**2) * (Q**2) * (gamma**2)
        
        
        term_comun = (numerador / denominador) * np.sinh(Q * L)**2
        
    else:
        
        Q = np.sqrt(C2 * (E - V0))
        
        if Q == 0: 
            term_comun = 0 
        else:
            
            numerador = ((K * gamma)**2 - Q**2)**2
            denominador = 4 * (K**2) * (Q**2) * (gamma**2)
            
            term_comun = (numerador / denominador) * np.sin(Q * L)**2
    
  
    T_val = 1 / (1 + term_comun)
    R_val = term_comun / (1 + term_comun)
    
    T_lista.append(T_val)
    R_lista.append(R_val)
    Suma_lista.append(T_val + R_val)


# E_n = (pi^2 * n^2) / (C2 * L^2) + V0
resonancias_E = []
resonancias_T = []

print("Resonancias Teóricas:")
for n in range(1, 5):
    E_res = (np.pi**2 * n**2) / (C2 * L**2) + V0
    
    if E_res <= max(Energias):
        resonancias_E.append(E_res)
        resonancias_T.append(1.0)
        print(f"n={n}: {E_res:.4f} eV")


color_c1 = '#E06C9F'  
color_c2 = '#5DADE2'  
color_c3 = '#A9DFBF'
color_c4 = '#F7DC6F'

plt.figure(figsize=(10, 7))

plt.plot(Energias, T_lista, label='Transmitancia ($T$)', color=color_c1, linewidth=2)
plt.plot(Energias, R_lista, label='Reflectancia ($R$)', color=color_c2, linewidth=2, linestyle='-')
plt.plot(Energias, Suma_lista, label='$T + R$ (Conservación)', color=color_c3, linewidth=1.5, alpha=0.7)


plt.plot(resonancias_E, resonancias_T, marker='o', linestyle='None', 
         label=r'Resonancias Teóricas ($E_\eta$)', color=color_c4, markersize=8, zorder=5)


plt.axvline(x=V0, color='gray', linestyle=':', label=f'Barrera $V_0$')


plt.title(f'Barrera Rectangular: (Unidades de Estado Sólido)\n$L={L}$ nm, $\\gamma={gamma:.2f}$', fontsize=14)
plt.xlabel('Energía (eV)', fontsize=12)
plt.ylabel('Coeficiente de Transmisión T', fontsize=12)
plt.legend(loc='center right')
plt.grid(True, alpha=0.3)
plt.ylim(-0.1, 1.2) 

plt.show()