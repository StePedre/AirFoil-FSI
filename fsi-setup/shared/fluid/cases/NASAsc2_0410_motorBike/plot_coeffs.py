import matplotlib.pyplot as plt
import pandas as pd

# Leggi il file saltando le righe di commento (#)
df = pd.read_csv("postProcessing/forceCoeffs1/0/coefficient_0.dat", 
                 delim_whitespace=True,
                 comment='#',
                 header=None)

# Dai nome alle colonne
df.columns = [
    "Time","Cd","Cd_f","Cd_r","Cl","Cl_f","Cl_r",
    "CmPitch","CmRoll","CmYaw",
    "Cs","Cs_f","Cs_r"
]

# Plot Cd
plt.figure()
plt.plot(df["Time"], df["Cd"])
plt.xlabel("Time")
plt.ylabel("Cd")
plt.title("Drag Coefficient")
plt.grid(True)

# Plot Cl
plt.figure()
plt.plot(df["Time"], df["Cl"])
plt.xlabel("Time")
plt.ylabel("Cl")
plt.title("Lift Coefficient")
plt.grid(True)

plt.show()
