import os
import re

def get_case_properties():
    # Percorsi dei file relativi alla cartella corrente
    current_dir = os.getcwd()
    snappy_path = os.path.join(current_dir, "system", "snappyHexMeshDict")
    decomp_path = os.path.join(current_dir, "system", "decomposeParDict")


    # ==========================================
    # 1. LETTURA MESH (Riga 18 di snappyHexMeshDict)
    # ==========================================
    mesh_info = "N/A (File mancante o errore lettura)"
    
    if os.path.exists(snappy_path):
        try:
            with open(snappy_path, 'r') as f:
                lines = f.readlines()
                # In Python gli indici partono da 0, quindi riga 18 = indice 17
                target_line_idx = 17 
                
                if len(lines) >= target_line_idx + 1:
                    raw_line = lines[target_line_idx].strip()
                    # Rimuovo i caratteri di commento "//" per pulire l'output
                    clean_text = raw_line.replace("//", "").strip()
                    mesh_info = clean_text if clean_text else "Riga vuota"
                else:
                    mesh_info = f"Errore: Il file ha meno di {target_line_idx+1} righe."
        except Exception as e:
            mesh_info = f"Errore lettura: {e}"
    else:
        mesh_info = "ERRORE: system/snappyHexMeshDict non trovato."

    # ==========================================
    # 2. LETTURA CPU (Riga 17 di decomposeParDict)
    # ==========================================
    cpu_info = "N/A (File mancante o errore lettura)"

    if os.path.exists(decomp_path):
        try:
            with open(decomp_path, 'r') as f:
                lines = f.readlines()
                # Riga 17 = indice 16
                target_line_idx = 16
                
                if len(lines) >= target_line_idx + 1:
                    raw_line = lines[target_line_idx].strip()
                    
                    # Cerchiamo di estrarre solo il numero usando Regex
                    # Cerca: numberOfSubdomains seguito da un numero
                    match = re.search(r'numberOfSubdomains\s+(\d+);', raw_line)
                    if match:
                        cpu_info = match.group(1) # Prende solo il numero (es. 16)
                    else:
                        # Se il formato è strano, stampa tutta la riga pulita
                        cpu_info = raw_line.replace(";", "")
                else:
                    cpu_info = f"Errore: Il file ha meno di {target_line_idx+1} righe."
        except Exception as e:
            cpu_info = f"Errore lettura: {e}"
    else:
        cpu_info = "ERRORE: system/decomposeParDict non trovato."

    # ==========================================
    # 3. STAMPA RISULTATI
    # ==========================================
    print(f"Mesh Type:  {mesh_info}")
    print(f"CPU Cores:  {cpu_info}") 

if __name__ == "__main__":
    get_case_properties()
