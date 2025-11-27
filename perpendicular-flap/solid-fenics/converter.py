import meshio
import numpy as np

# 1. Leggi il file originale di Gmsh
print("Lettura file .msh...")
mesh = meshio.read("NASAsc2-0410.msh")

# 2. Estrai i TETRAEDRI (Volume) - Linearizza se necessario
print("Estrazione Volume (Tetra)...")
tetra_cells = []
for cell_block in mesh.cells:
    if cell_block.type == "tetra10": # Quadratic -> Linear
        tetra_cells.append(("tetra", cell_block.data[:, :4]))
    elif cell_block.type == "tetra":
        tetra_cells.append(("tetra", cell_block.data))

# Crea mesh volume pulita
volume_mesh = meshio.Mesh(points=mesh.points, cells=tetra_cells)
# Scrivi in XML (FEniCS nativo)
meshio.write("mesh.xml", volume_mesh)

# 3. Estrai i TRIANGOLI (Boundary) e i loro TAG
print("Estrazione Boundary (Triangle)...")
triangle_cells = []
triangle_data = []

# Cerca i dati fisici (tags)
key = "gmsh:physical"
# Fallback se la chiave è diversa
if key not in mesh.cell_data:
    if len(mesh.cell_data) > 0:
        key = list(mesh.cell_data.keys())[0]

for cell_block, data in zip(mesh.cells, mesh.cell_data[key]):
    if cell_block.type == "triangle6": # Quadratic -> Linear
        triangle_cells.append(("triangle", cell_block.data[:, :3]))
        triangle_data.append(data)
    elif cell_block.type == "triangle":
        triangle_cells.append(("triangle", cell_block.data))
        triangle_data.append(data)

# Crea mesh boundary temporanea
# Nota: "name_to_read" è il nome interno che daremo ai dati
boundary_mesh = meshio.Mesh(
    points=mesh.points, 
    cells=triangle_cells, 
    cell_data={"name_to_read": triangle_data}
)

# Scrivi in XML (Questo crea boundaries.xml E boundaries_name_to_read.xml)
# boundaries_name_to_read.xml è quello che contiene la MeshFunction che ci serve!
meshio.write("boundaries.xml", boundary_mesh)

print("------------------------------------------------")
print("✅ Fatto! File generati:")
print("   - mesh.xml (Volume)")
print("   - boundaries_name_to_read.xml (Usa questo in solid.py!)")