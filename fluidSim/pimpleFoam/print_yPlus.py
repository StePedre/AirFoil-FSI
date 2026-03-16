import sys
import pandas as pd
import os

def main():
    # 1. Controllo argomenti da riga di comando
    if len(sys.argv) < 2:
        print("Errore: Specifica il percorso del file yPlus.")
        print("Uso: python print_latest_yplus.py /percorso/al/file/yPlus.dat")
        sys.exit(1)

    file_path = sys.argv[1]

    if not os.path.exists(file_path):
        print(f"Errore: Il file '{file_path}' non esiste.")
        sys.exit(1)

    # 2. Lettura del file
    # I file OpenFOAM usano '#' per i commenti e spazi multipli come separatori.
    # Definiamo i nomi delle colonne manualmente perché l'header nel file inizia con #.
    column_names = ['Time', 'patch', 'min', 'max', 'average']

    try:
        df = pd.read_csv(
            file_path,
            sep=r'\s+',       # Espressione regolare per separatori spazi multipli
            comment='#',      # Ignora le righe che iniziano con #
            names=column_names,
            header=None       # Non c'è header (perché è commentato)
        )
    except Exception as e:
        print(f"Errore durante la lettura del file: {e}")
        sys.exit(1)

    if df.empty:
        print("Il file è vuoto o non contiene dati validi.")
        sys.exit(0)

    # 3. Identificazione dell'ultimo timestep
    # Trova il valore massimo nella colonna 'Time'
    latest_time = df['Time'].max()

    # Filtra il DataFrame mantenendo solo le righe corrispondenti all'ultimo tempo
    latest_data = df[df['Time'] == latest_time]

    # 4. Stampa dei risultati
    print("-" * 54)
    print(f" ANALISI Y+ PER L'ULTIMO TIMESTEP (Time = {latest_time} s)")
    print("-" * 54)
    # Intestazione tabella
    print(f" {'PATCH':<25} | {'AVG y+':<15}")
    print("-" * 54)

    # Itera su ogni patch trovata nell'ultimo timestep
    for _, row in latest_data.iterrows():
        patch_name = row['patch']
        avg_y = row['average']

        # Stampa formattata
        print(f" {patch_name:<25} | {avg_y:<15.4f}")

    print("-" * 54)

if __name__ == "__main__":
    main()
