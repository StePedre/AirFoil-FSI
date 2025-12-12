import os
import shutil
import argparse
import sys

def setup_simulation(angle, mesh_type):
    # ==========================================
    # 1. CONFIGURAZIONE PERCORSI RELATIVI
    # ==========================================
    
    # Cartella corrente (dove gira lo script, es: .../NASAsc2-0410_pimpleFoam)
    current_case_path = os.getcwd()
    
    # Percorsi delle librerie esterne
    meshes_library_path = os.path.join(current_case_path, "../../meshes")
    stls_library_path = os.path.join(current_case_path, "../../../cadFiles")

    # ==========================================
    # 2. LOGICA SELEZIONE FILE STL
    # ==========================================
    
    # Mappa i nomi dei file STL in base all'angolo
    stl_map = {
        0:  "NASAsc2-0410_singleLine.stl",
        5:  "NASAsc2-0410_singleLine_5deg.stl",
        10: "NASAsc2-0410_singleLine_10deg.stl"
    }
    
    source_stl_filename = stl_map[angle]
    
    # Percorso completo sorgente STL
    src_stl_path = os.path.join(stls_library_path, source_stl_filename)
    
    # Percorso destinazione STL (Rinominiamo sempre in un nome standard per OpenFOAM)
    dst_stl_path = os.path.join(current_case_path, "constant", "triSurface", "NASAsc2-0410.stl")

    # ==========================================
    # 3. LOGICA SELEZIONE FILE MESH DICT
    # ==========================================
    
    # Nome cartella preset (es: 0coarse)
    preset_folder_name = f"{angle}{mesh_type}"
    
    # Percorso completo sorgente Dict
    src_dict_path = os.path.join(meshes_library_path, preset_folder_name, "snappyHexMeshDict")
    
    # Percorso destinazione Dict
    dst_dict_path = os.path.join(current_case_path, "system", "snappyHexMeshDict")

    # ==========================================
    # 4. CONTROLLI DI SICUREZZA
    # ==========================================
    
    print(f"--- Setup Caso: Angolo {angle}° | Mesh {mesh_type} ---")
    
    # Verifica esistenza file STL sorgente
    if not os.path.exists(src_stl_path):
        print(f"\n[ERRORE] File STL non trovato: {src_stl_path}")
        print(f"Verifica la cartella '../stlFiles/'")
        sys.exit(1)

    # Verifica esistenza file Dict sorgente
    if not os.path.exists(src_dict_path):
        print(f"\n[ERRORE] snappyHexMeshDict non trovato: {src_dict_path}")
        print(f"Verifica la cartella '../meshes/{preset_folder_name}/'")
        sys.exit(1)

    # Crea le cartelle di destinazione se mancano
    os.makedirs(os.path.dirname(dst_stl_path), exist_ok=True)
    os.makedirs(os.path.dirname(dst_dict_path), exist_ok=True)

    # ==========================================
    # 5. ESECUZIONE COPIA
    # ==========================================
    
    try:
        # Copia STL
        shutil.copy2(src_stl_path, dst_stl_path)
        print(f"Copiato STL:  ../stlFiles/{source_stl_filename}\n             -> constant/triSurface/NASAsc2-0410.stl")

        # Copia Dict
        shutil.copy2(src_dict_path, dst_dict_path)
        print(f"Copiato Dict: ../meshes/{preset_folder_name}/snappyHexMeshDict\n             -> system/snappyHexMeshDict")
        
        print("\n--> Setup completato.")

    except Exception as e:
        print(f"\n[CRITICO] Errore di copia file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Configura il caso copiando STL da stlFiles e Dict da meshes.")
    
    # Flag Angolo
    parser.add_argument(
        "-a", "--angle", 
        type=int, 
        choices=[0, 5, 10], 
        required=True,
        help="Angolo di attacco (0, 5, 10)"
    )
    
    # Flag Tipo Mesh
    parser.add_argument(
        "-m", "--mesh", 
        type=str, 
        choices=["coarse", "refined"], 
        required=True,
        help="Tipo di mesh (coarse, refined)"
    )

    args = parser.parse_args()
    
    setup_simulation(args.angle, args.mesh)
