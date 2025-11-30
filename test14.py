import numpy as np
import matplotlib.pyplot as plt

nm_to_bohr = 18.8973        # 1 nm = ~18.9 bohr
eV_to_Hartree = 0.0367493   # 1 eV = ~0.037 Ha

L_nm = 10.0      #nm    
V0_eV = 0.23        # eV
m1_eff = 0.067       
m2_eff = 0.092       

L_au = L_nm * nm_to_bohr      
V0_au = V0_eV * eV_to_Hartree 

Energias_eV = np.linspace(0.001, 0.6, 1000)
Energias_Ha = Energias_eV * eV_to_Hartree


gamma = m2_eff / m1_eff


C1 = 2 * m1_eff 
C2 = 2 * m2_eff


T_lista = []
R_lista = []
Suma_lista = []

for E in Energias_Ha: 
    
    # K y Q resultan en unidades de [1/Bohr]
    K = np.sqrt(C1 * E)
    
    term_comun = 0.0
    
    if E < V0_au:
        Q = np.sqrt(C2 * (V0_au - E))
        
        numerador = ((K * gamma)**2 + Q**2)**2
        denominador = 4 * (K**2) * (Q**2) * (gamma**2)
        
        
        term_comun = (numerador / denominador) * np.sinh(Q * L_au)**2
        
    else:
        Q = np.sqrt(C2 * (E - V0_au))
        
        if Q == 0: 
            term_comun = 0 
        else:
            numerador = ((K * gamma)**2 - Q**2)**2
            denominador = 4 * (K**2) * (Q**2) * (gamma**2)
            
            term_comun = (numerador / denominador) * np.sin(Q * L_au)**2
    
    T_val = 1 / (1 + term_comun)
    R_val = term_comun / (1 + term_comun)
    
    T_lista.append(T_val)
    R_lista.append(R_val)
    Suma_lista.append(T_val + R_val)



resonancias_E_eV = [] 
resonancias_T = []

print("Resonancias Teóricas :")
for n in range(1, 10):
    
    E_res_au = (np.pi**2 * n**2) / (C2 * L_au**2) + V0_au
    
    
    E_res_eV = E_res_au / eV_to_Hartree
    
    if E_res_eV <= max(Energias_eV):
        resonancias_E_eV.append(E_res_eV)
        resonancias_T.append(1.0)
        print(f"n={n}: {E_res_eV:.4f} eV")

color_c1 = "#1ABD0F"  
color_c2 = "#A800BE"  
color_c3 = "#AD1919"
color_c4 = "#938CBE"

plt.figure(figsize=(10, 7))


plt.plot(Energias_eV, T_lista, label='Transmitancia ($T$)', color=color_c1, linewidth=2)
plt.plot(Energias_eV, R_lista, label='Reflectancia ($R$)', color=color_c2, linewidth=2, linestyle='-')
plt.plot(Energias_eV, Suma_lista, label='$T + R$', color=color_c3, linewidth=1.5, alpha=0.7)

plt.plot(resonancias_E_eV, resonancias_T, marker='o', linestyle='None', 
         label=r'Resonancias ($E_n$)', color=color_c4, markersize=8, zorder=5)

plt.axvline(x=V0_eV, color='gray', linestyle=':', label=f'Barrera $V_0$')


plt.title(f'Barrera Rectangular (Unidades Atómicas)\n$L={L_nm}$ nm, $\\gamma={gamma:.2f}$', fontsize=14)
plt.xlabel('Energía (Ha)', fontsize=12)
plt.ylabel('Coeficiente de Transmisión', fontsize=12)
plt.legend(loc='center right')
plt.grid(True, alpha=0.3)
plt.ylim(-0.1, 1.2) 

plt.show()