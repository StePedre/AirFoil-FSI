import sys
import numpy as np

# --- Configurazione Utente ---
# Numero di righe (passi temporali) da IGNORARE all'inizio.
# Questi primi passi rappresentano la fase di transitorio iniziale (start-up).
# Basandosi sui tuoi dati di esempio (molto brevi), impostiamo un valore conservativo.
ROWS_TO_SKIP = 50 

# Colonna di Cd (indice 1) e Cl (indice 4) nel tuo file.
# [Time, Cd, Cd(f), Cd(r), Cl, Cl(f), ...]
CD_COL = 1
CL_COL = 4

def calculate_mean_coeffs():
    """
    Reads forceCoeffs data, skips the initial transient phase,
    and calculates the temporal average of Cd and Cl.
    """
    if len(sys.argv) < 2:
        print("Error: Please provide the file path as a command-line argument.")
        print("Usage: python calculate_average_coeffs.py /path/to/forceCoeffs.dat")
        sys.exit(1)

    file_path = sys.argv[1]
    
    # OpenFOAM forceCoeffs files use arbitrary whitespace as separator.
    # The header starts with '#' and needs to be skipped.
    try:
        # Read all data, skip the first row (the header starting with '#')
        # We use header=None and skip initial rows for robust reading.
        df = pd.read_csv(
            file_path,
            sep=r'\s+',  # Regex for one or more whitespace characters
            comment='#', # Ignore commented lines (including the header)
            header=None,
            skipinitialspace=True
        )
    except FileNotFoundError:
        print(f"Error: File not found at the specified path: {file_path}")
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred while reading the file: {e}")
        sys.exit(1)

    # Assicurati che il DataFrame non sia vuoto e pulisci i NaN (a causa degli spazi multipli)
    df = df.dropna(axis=1, how='all')
    
    # 1. Saltare la fase di transitorio iniziale
    if len(df) <= ROWS_TO_SKIP:
        print(f"Warning: Only {len(df)} time steps found. Cannot skip {ROWS_TO_SKIP} steps for averaging.")
        start_row = 0
    else:
        start_row = ROWS_TO_SKIP
        
    df_stable = df.iloc[start_row:]
    
    if df_stable.empty:
        print("Error: DataFrame is empty after filtering the transient phase.")
        sys.exit(1)
        
    # 2. Estrai le colonne di Cd e Cl (basato sullo schema di numerazione del file, indice 0 = Time)
    # OpenFOAM files often list Time as the first column, so Cd is index 1, Cl is index 4.
    try:
        Cd_values = df_stable.iloc[:, CD_COL].astype(float)
        Cl_values = df_stable.iloc[:, CL_COL].astype(float)
    except IndexError:
        print("Error: Column index out of bounds. Check if the file structure has changed.")
        print(f"Expected at least {CL_COL + 1} columns.")
        sys.exit(1)

    # 3. Calcola la media temporale
    avg_Cd = np.mean(Cd_values)
    avg_Cl = np.mean(Cl_values)
    
    # 4. Stampa i risultati
    print("--- Temporal Averaging Results ---")
    print(f"Total time steps analyzed: {len(df)} (from {df.iloc[0, 0]:.4f}s to {df.iloc[-1, 0]:.4f}s)")
    print(f"Skipped initial {start_row} steps (Transient Phase).")
    print(f"Averaging over {len(df_stable)} stable time steps.")
    print("-" * 35)
    print(f"Average Drag Coefficient (Cd): {avg_Cd:.6e}")
    print(f"Average Lift Coefficient (Cl): {avg_Cl:.6e}")
    print("-" * 35)

if __name__ == "__main__":
    calculate_mean_coeffs()
