import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq

# --- CONFIGURAZIONE ---
filename = "precice-Solid-watchpoint-Flap-Tip.log"
dt = 0.0025  # Assicurati che corrisponda al tuo Dt di simulazione

# --- CARICAMENTO DATI ---
# Saltiamo la prima riga di intestazione
data = np.loadtxt(filename, skiprows=1)

time = data[:, 0]
# Coordinate iniziali (colonne 1, 2, 3)
# Spostamento (colonne 4, 5, 6) -> 4=X, 5=Y (Portanza), 6=Z
disp_y = data[:, 5] * 1000  # Convertiamo in mm
# Forze (colonne 7, 8, 9) -> 8=Y (Forza verticale)
force_y = data[:, 8]

# --- 1. PLOT SPOSTAMENTO NEL TEMPO ---
plt.figure(figsize=(12, 8))

plt.subplot(3, 1, 1)
plt.plot(time, disp_y, label='Spostamento Verticale (Y)', color='blue', linewidth=1.5)
plt.title('Analisi FSI Ala: Spostamento Punta')
plt.ylabel('Spostamento [mm]')
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend()

# --- 2. PLOT FORZA VERTICALE (LIFT LOCALE) ---
plt.subplot(3, 1, 2)
plt.plot(time, force_y, label='Forza Fluida (Y)', color='red', linewidth=1.5)
plt.ylabel('Forza [N]')
plt.xlabel('Tempo [s]')
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend()

# --- 3. ANALISI IN FREQUENZA (FFT) ---
# Calcoliamo la FFT per trovare la frequenza propria accoppiata
N = len(time)
yf = fft(disp_y - np.mean(disp_y)) # Sottraiamo la media per eliminare la componente statica
xf = fftfreq(N, dt)[:N//2]

plt.subplot(3, 1, 3)
plt.plot(xf, 2.0/N * np.abs(yf[0:N//2]), color='green')
plt.title('Spettro di Frequenza (FFT)')
plt.xlabel('Frequenza [Hz]')
plt.ylabel('Ampiezza')
plt.xlim(0, 50) # Limitiamo a 50Hz per vedere meglio i primi modi
plt.grid(True, linestyle='--', alpha=0.7)

plt.tight_layout()
plt.savefig('watchpoint_analisys.png')
plt.show()

# --- STATISTICHE ---
print(f"Spostamento Massimo: {np.max(disp_y):.4f} mm")
print(f"Spostamento Medio (Equilibrio): {np.mean(disp_y):.4f} mm")
peak_freq = xf[np.argmax(2.0/N * np.abs(yf[0:N//2]))]
print(f"Frequenza di Picco rilevata: {peak_freq:.2f} Hz")