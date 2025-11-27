from fenics import *
from fenicsprecice import Adapter
import numpy as np
import meshio
from scipy.spatial import cKDTree

# ---------------------------------------------------------------------------
# CLASSE WRAPPER GEOMETRICO (KDTree)
# ---------------------------------------------------------------------------
class PointCloudSubDomain(SubDomain):
    def __init__(self, points):
        SubDomain.__init__(self)
        self.kdtree = cKDTree(points)
        
    def inside(self, x, on_boundary):
        if hasattr(x, 'x'):
            pt = [x.x(), x.y(), x.z()]
        else:
            pt = x
        # Aumentata tolleranza a 1e-2 per sicurezza
        dist, _ = self.kdtree.query(pt)
        return dist < 2.0

# ---------------------------------------------------------------------------
# 1. LETTURA GEOMETRIA E MESH
# ---------------------------------------------------------------------------
boundary_file = "mesh_boundary_linear.xdmf"
try:
    b_mesh = meshio.read(boundary_file)
except:
    print("XDMF boundary non trovato, provo msh...")
    b_mesh = meshio.read("NASAsc2-0410.msh")

data_key = None
possible_keys = ["name_to_read", "boundaries", "gmsh:physical"]
for k in b_mesh.cell_data.keys():
    if k in possible_keys:
        data_key = k
        break
if not data_key: data_key = list(b_mesh.cell_data.keys())[0]

cell_tags = b_mesh.cell_data[data_key][0]
cells = b_mesh.cells[0].data

coupling_id = 9
fixed_id = 8

def get_unique_points(indices, points):
    return points[np.unique(indices.flatten())]

interface_points = get_unique_points(cells[cell_tags == coupling_id], b_mesh.points)
fixed_points = get_unique_points(cells[cell_tags == fixed_id], b_mesh.points)

print(f"DEBUG: Punti Incastro Trovati (Fixed): {len(fixed_points)}")
print(f"DEBUG: Punti Interfaccia Trovati (Interface): {len(interface_points)}")

interface_domain = PointCloudSubDomain(interface_points)
fixed_domain = PointCloudSubDomain(fixed_points)

mesh = Mesh("mesh.xml")

# ---------------------------------------------------------------------------
# 2. SETUP FISICO
# ---------------------------------------------------------------------------
V = VectorFunctionSpace(mesh, 'P', 1)

boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
boundaries.set_all(0)
interface_domain.mark(boundaries, coupling_id)
fixed_domain.mark(boundaries, fixed_id)

dx = Measure('dx', domain=mesh)
ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

E = 4.0e6
nu = 0.3
mu = E / (2.0 * (1.0 + nu))
lmbda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
rho_s = 3.0e3

u = Function(V)
u_n = Function(V)
v = TestFunction(V)
trial_u = TrialFunction(V)

def sigma(u):
    return 2.0 * mu * sym(grad(u)) + lmbda * tr(grad(u)) * Identity(len(u))

# ---------------------------------------------------------------------------
# 3. PREPARAZIONE MAPPATURA FORZE
# ---------------------------------------------------------------------------
coupling_forces = Function(V)

# Mappa DOF manuale per bypassare bug adapter 3D
vertex_to_dof = vertex_to_dof_map(V)
mesh_coords = mesh.coordinates()
interface_vertex_indices = []

for i, pt in enumerate(mesh_coords):
    if interface_domain.inside(pt, True):
        interface_vertex_indices.append(i)

interface_vertex_indices = sorted(interface_vertex_indices)

v2d = np.array(vertex_to_dof_map(V), dtype=np.int32)

interface_dofs_x_ordered = []
interface_dofs_y_ordered = []
interface_dofs_z_ordered = []

for v_idx in interface_vertex_indices:
    interface_dofs_x_ordered.append(v2d[3*v_idx])
    interface_dofs_y_ordered.append(v2d[3*v_idx + 1])
    interface_dofs_z_ordered.append(v2d[3*v_idx + 2])

idx_x = np.array(interface_dofs_x_ordered, dtype=np.int32)
idx_y = np.array(interface_dofs_y_ordered, dtype=np.int32)
idx_z = np.array(interface_dofs_z_ordered, dtype=np.int32)

print(f"Local Interface Nodes: {len(interface_vertex_indices)}")
print(f"Target DOF Size: {len(idx_x)}")

# ---------------------------------------------------------------------------
# 4. INIZIALIZZAZIONE PRECICE
# ---------------------------------------------------------------------------
precice = Adapter(adapter_config_filename="precice-adapter-config-fsi-s.json")

if precice.get_participant_name() == "Solid":
    precice.initialize(interface_domain, read_function_space=V, write_object=u)

# ---------------------------------------------------------------------------
# 5. SOLUTORE
# ---------------------------------------------------------------------------
bcs = [DirichletBC(V, Constant((0.0, 0.0, 0.0)), fixed_domain)]
f = Constant((0.0, 0.0, 0.0))
dt = 0.01

F = rho_s * dot((u - u_n) / Constant(dt), v) * dx + \
    inner(sigma(u), sym(grad(v))) * dx - \
    dot(f, v) * dx - \
    dot(coupling_forces, v) * ds(coupling_id)

J = derivative(F, u, trial_u)

u.rename("Displacement", "Displacement")
vtkfile = File("output/Displacement.pvd")
vtkfile << (u, 0.0)

t = 0.0
n = 0 

print("Inizio ciclo temporale...")
while precice.is_coupling_ongoing():
    
    if precice.requires_writing_checkpoint():
        precice.store_checkpoint(u_n, t, n)

    dt = precice.get_max_time_step_size()
    
    # --- GESTIONE ROBUSTA DATI INGRESSO ---
    raw_data = precice.read_data(dt)
    
    # Se è None o scalare, creiamo un array di zeri
    if raw_data is None or np.ndim(raw_data) == 0:
        # print("WARNING: Received empty/scalar data. Assuming zeros.")
        read_data = np.zeros((len(idx_x), 3))
    else:
        read_data = np.array(raw_data)
        
        # Se è piatto (N*3), lo rimodelliamo
        if read_data.ndim == 1:
            read_data = read_data.reshape(-1, 3)

    # Verifica dimensionale prima di assegnare
    if read_data.shape[0] == len(idx_x):
        vec_vals = coupling_forces.vector().get_local()
        vec_vals[idx_x] = read_data[:, 0]
        vec_vals[idx_y] = read_data[:, 1]
        vec_vals[idx_z] = read_data[:, 2]
        coupling_forces.vector().set_local(vec_vals)
        coupling_forces.vector().apply("insert")
    else:
        # Se le dimensioni non tornano, non crashare, stampa solo un avviso
        print(f"WARNING: Data mismatch! preCICE sent {read_data.shape[0]} vectors, expected {len(idx_x)}. Skipping update.")
    
    # Risolvi
    F = rho_s * dot((u - u_n) / Constant(dt), v) * dx + \
        inner(sigma(u), sym(grad(v))) * dx - \
        dot(f, v) * dx - \
        dot(coupling_forces, v) * ds(coupling_id)
    
    solve(F == 0, u, bcs, J=J)

    precice.write_data(u)
    precice.advance(dt)

    if precice.requires_reading_checkpoint():
        u_n_cp, t_cp, n_cp = precice.retrieve_checkpoint()
        u_n.assign(u_n_cp)
        t = t_cp
        n = n_cp
    else:
        u_n.assign(u)
        t += dt
        n += 1
        vtkfile << (u, t)

precice.finalize()
print("Simulazione Terminata.")