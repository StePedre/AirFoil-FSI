import pandas as pd
import matplotlib.pyplot as plt
import os

# ==========================================
# PARAMETRI FISICI (Modifica qui)
# ==========================================
rho = 1.225        # Densità dell'aria [kg/m^3]
v = 100.0           # Velocità indisturbata [m/s]
S = 3            # Superficie alare [m^2]

file_path = "postProcessing/forces/0/coefficient.dat"
output_coeffs = "plot_coefficienti.png"
output_forces = "plot_forze.png"

# ==========================================
# CARICAMENTO E PREPARAZIONE DATI
# ==========================================
if not os.path.exists(file_path):
    print(f"Errore: File {file_path} non trovato!")
else:
    # Lettura dati
    df = pd.read_csv(file_path, sep=r'\s+', comment='#', header=None)
    df.columns = ['Time', 'Cd', 'Cd_f', 'Cd_r', 'Cl', 'Cl_f', 'Cl_r', 'CmPitch', 'CmRoll', 'CmYaw', 'Cs', 'Cs_f', 'Cs_r']

    # Calcolo Forze
    q = 0.5 * rho * (v**2)
    df['Lift_N'] = q * S * df['Cl']
    df['Drag_N'] = q * S * df['Cd']

    # ------------------------------------------
    # GRAFICO 1: COEFFICIENTI (Cl, Cd)
    # ------------------------------------------
    plt.figure(figsize=(10, 6))
    plt.plot(df['Time'], df['Cl'], label='Cl (Portanza)', color='blue')
    plt.plot(df['Time'], df['Cd'], label='Cd (Resistenza)', color='red')
    plt.xlabel('Tempo/Iterazioni')
    plt.ylabel('Coefficienti [-]')
    plt.title('Andamento Coefficienti Aerodinamici')
    plt.legend()
    plt.grid(True, linestyle='--')
    
    plt.savefig(output_coeffs, dpi=300)
    plt.close() # Chiude la figura senza mostrarla
    print(f"Grafico coefficienti salvato: {output_coeffs}")

    # ------------------------------------------
    # GRAFICO 2: FORZE IN NEWTON (L, D)
    # ------------------------------------------
    plt.figure(figsize=(10, 6))
    plt.plot(df['Time'], df['Lift_N'], label='Lift (Portanza)', color='darkblue', linestyle='--')
    plt.plot(df['Time'], df['Drag_N'], label='Drag (Resistenza)', color='darkred', linestyle='--')
    plt.xlabel('Tempo/Iterazioni')
    plt.ylabel('Forza [N]')
    plt.title('Andamento Forze Aerodinamiche')
    plt.legend()
    plt.grid(True, linestyle='--')
    
    plt.savefig(output_forces, dpi=300)
    plt.close() # Chiude la figura senza mostrarla
    print(f"Grafico forze salvato: {output_forces}")

    # ------------------------------------------
    # STAMPA RIASSUNTO IN CONSOLE
    # ------------------------------------------
    last_100 = df.iloc[-100:]
    print("\n--- VALORI MEDI FINALI ---")
    print(f"Cl: {last_100['Cl'].mean():.4f} | Cd: {last_100['Cd'].mean():.4f}")
    print(f"Portanza: {last_100['Lift_N'].mean():.2f} N")
    print(f"Resistenza: {last_100['Drag_N'].mean():.2f} N")