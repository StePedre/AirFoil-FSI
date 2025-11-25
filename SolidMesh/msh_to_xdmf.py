import meshio
import numpy as np

# --- CONFIGURAZIONE ---
input_file = "NASAsc2-0410.msh"
domain_output = "mesh_domain.xdmf"
boundary_output = "mesh_boundary.xdmf"

print(f"Leggo {input_file}...")
try:
    msh = meshio.read(input_file)
except ValueError as e:
    print("\nERRORE CRITICO DI LETTURA:")
    print(e)
    print("SOLUZIONE: In Gmsh, esporta come 'Version 2.2 ASCII' e togli la spunta a 'Save all elements'.\n")
    exit()

# Funzione per estrarre celle e dati per un tipo specifico (es. 'tetra' o 'triangle')
def create_mesh_for_type(mesh_obj, cell_type_str):
    cells = []
    cell_data = []
    
    # Itera su tutti i blocchi trovati da meshio
    for i, cell_block in enumerate(mesh_obj.cells):
        # Controlla se il blocco è del tipo che cerchiamo (es. "tetra")
        # Gestisce anche varianti come "tetra10" (ordine 2)
        if cell_type_str in cell_block.type:
            cells.append(cell_block.data)
            # Cerca i dati fisici corrispondenti (se esistono)
            if "gmsh:physical" in mesh_obj.cell_data_dict:
                cell_data.append(mesh_obj.cell_data_dict["gmsh:physical"][cell_block.type])
            elif "gmsh:physical" in mesh_obj.cell_data: # Supporto legacy meshio
                cell_data.append(mesh_obj.cell_data["gmsh:physical"][i])

    if not cells:
        return None

    # Unisce tutti i blocchi trovati in un unico blocco massiccio
    cells_concatenated = np.concatenate(cells)
    
    mesh_args = {"points": mesh_obj.points, "cells": {cell_type_str: cells_concatenated}}
    
    if cell_data:
        data_concatenated = np.concatenate(cell_data)
        # "name_to_read" è un nome fittizio che FEniCS userà per leggere i tag
        mesh_args["cell_data"] = {"name_to_read": [data_concatenated]}
    
    return meshio.Mesh(**mesh_args)

# --- 1. VOLUME (Cerca tetraedri) ---
# Proviamo a cercare tetraedri (tetra) o tetraedri quadratici (tetra10)
element_target = "tetra"
# Se la mesh è ordine 2, meshio vede "tetra10", altrimenti "tetra"
# Questo trucco trova il nome esatto presente nel file
for cell in msh.cells:
    if "tetra" in cell.type:
        element_target = cell.type
        break

print(f"Estraggo il volume (tipo elementi: {element_target})...")
mesh_vol = create_mesh_for_type(msh, element_target)

if mesh_vol:
    meshio.write(domain_output, mesh_vol)
    print(f"-> Volume salvato in {domain_output}")
else:
    print("ERRORE: Nessun tetraedro trovato! Hai creato il Physical Volume in Gmsh?")

# --- 2. SUPERFICI (Cerca triangoli) ---
surf_target = "triangle"
for cell in msh.cells:
    if "triangle" in cell.type:
        surf_target = cell.type
        break

print(f"Estraggo le superfici (tipo elementi: {surf_target})...")
mesh_surf = create_mesh_for_type(msh, surf_target)

if mesh_surf:
    meshio.write(boundary_output, mesh_surf)
    print(f"-> Superfici salvate in {boundary_output}")   # CORRETTO
else:
    print("ATTENZIONE: Nessun triangolo trovato per le boundary.")

print("\nFinito!")