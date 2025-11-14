import matplotlib.pyplot as plt
import pandas as pd

file_path = "postProcessing/forceCoeffs1/0/coefficient_0.dat"

df = pd.read_csv(file_path, sep=r'\s+', comment='#', header=None)

df.columns = ['time', 'Cd', 'Cd(f)', 'Cd(r)',
              'Cl', 'Cl(f)', 'Cl(r)',
              'CmPitch', 'CmRoll', 'CmYaw',
              'Cs', 'Cs(f)', 'Cs(r)'][:df.shape[1]]

Cd = df['Cd'].to_numpy()
Cl = df['Cl'].to_numpy()

# Plot Cl vs Cd
plt.figure(figsize=(8, 6))
plt.plot(Cd, Cl, 'o-', markersize=3, color='blue')
plt.xlabel('Drag Coefficient (Cd)')
plt.ylabel('Lift Coefficient (Cl)')
plt.title('Lift vs Drag Coefficients')
plt.grid(True)
plt.tight_layout()
plt.show()

# Plot Cl vs time
plt.figure(figsize=(10, 6))
plt.plot(df['time'].to_numpy(), df['Cl'].to_numpy(), color='red', linewidth=2) 
plt.xlabel('Time')
plt.ylabel('Lift Coefficient (Cl)')
plt.title('Lift Coefficient (Cl) vs Time')
plt.grid(True)
plt.tight_layout()
plt.show()