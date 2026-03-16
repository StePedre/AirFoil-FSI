# Import required libs
from fenics import Constant, Function, AutoSubDomain, VectorFunctionSpace, interpolate, TrialFunction, TestFunction, Expression, DirichletBC, Identity, inner, dx, sym, grad, div, lhs, rhs, File, solve, assemble_system, Mesh, XDMFFile, dot, PETScKrylovSolver
import numpy as np
from fenicsprecice import Adapter
import time


# 1. Importa la Mesh (dominio 3D)
mesh = Mesh()
with XDMFFile("mesh/NASA_10deg_small.xdmf") as infile:
    infile.read(mesh)

n_vertices = mesh.num_vertices()
n_cells = mesh.num_cells()
print(f"\n[SOLID] Mesh caricata: {n_vertices} vertici, {n_cells} celle.")

# --- PROPRIETÀ DEL MATERIALE (Alluminio 7075) ---
dim = 3             # [adimensionale]
E = 71.7e9          # [Pa] o [N/m^2] - Modulo di Young
nu = 0.33           # [adimensionale] - Coefficiente di Poisson
rho = 2810.0        # [kg/m^3] - Densità dell'alluminio

g = 9.81            # [m/s^2] - Accelerazione di gravità
f_vol = Constant((0.0, -rho * g, 0.0)) # [N/m^3] - Forza di volume (peso)

# Costanti di Lamé
mu = Constant(E / (2.0 * (1.0 + nu)))
lambda_ = Constant(E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu)))

# create Function Space
V = VectorFunctionSpace(mesh, 'P', 1)

# Trial and Test Functions
du = TrialFunction(V)
v = TestFunction(V)

u_np1 = Function(V)

# function known from previous timestep
u_n = Function(V)
v_n = Function(V)
a_n = Function(V)

f_N_function = interpolate(Expression(("1", "0", "0"), degree=1), V)
u_function = interpolate(Expression(("0", "0", "0"), degree=1), V)


TOL = 1E-9 # # [m] tolleranza incastro

# Trova la Z minima della mesh
z_coordinates = mesh.coordinates()[:, 2] # Colonna Z
Z_MIN = np.min(z_coordinates)
Z_MAX = np.max(z_coordinates)

def clamped_boundary(x, on_boundary):
    return on_boundary and abs(x[2] - Z_MIN) < TOL

def neumann_boundary(x, on_boundary):
    return on_boundary and abs(x[2] - Z_MIN) >= TOL

fixed_boundary = AutoSubDomain(clamped_boundary)   
coupling_boundary = AutoSubDomain(neumann_boundary)

# --- Re-inizializzazione di BC e preCICE ---

# La BC di Dirichlet viene ora applicata tramite il Sottodominio geometrico.
bc = DirichletBC(V, Constant((0.0, 0.0, 0.0)), fixed_boundary)

print("[SOLID] Inizializzazione adattatore preCICE...")
precice = Adapter(adapter_config_filename="precice-adapter-config-fsi-s.json")

# Initialize the coupling interface
# preCICE userà coupling_boundary e fixed_boundary definiti geometricamente.
precice.initialize(coupling_boundary, read_function_space=V, write_object=V, fixed_boundary=fixed_boundary)
print("[SOLID] preCICE inizializzato correttamente.")

precice_dt = precice.get_max_time_step_size()
fenics_dt = precice_dt  # if fenics_dt == precice_dt, no subcycling is applied
# n_substeps = 5  # number of substeps per window
# fenics_dt = precice_dt / n_substeps  # if fenics_dt < precice_dt, subcycling is applied
dt = Constant(np.min([precice_dt, fenics_dt]))

"""
Check requirements for alpha_m and alpha_f from
    Chung, J., and Hulbert, G. M. (June 1, 1993). "A Time Integration Algorithm for Structural Dynamics With Improved Numerical Dissipation:
    The Generalized-α Method." ASME. J. Appl. Mech. June 1993; 60(2): 371–375. https://doi.org/10.1115/1.2900803
"""

# alpha method parameters
alpha_m = Constant(0.2)
alpha_f = Constant(0.4)

assert (float(alpha_m) <= float(alpha_f))
assert (float(alpha_f) <= 0.5)

gamma = Constant(0.5 + alpha_f - alpha_m)
beta = Constant((gamma + 0.5)**2 / 4.)


# Define strain
def epsilon(u):
    return 0.5 * (grad(u) + grad(u).T)


# Define Stress tensor
def sigma(u):
    return lambda_ * div(u) * Identity(dim) + 2 * mu * epsilon(u)


# Define Mass form
def m(u, v):
    return rho * inner(u, v) * dx


# Elastic stiffness form
def k(u, v):
    return inner(sigma(u), sym(grad(v))) * dx


# # External Work
# def Wext(u_):
#     return dot(u_, p) * ds
def W_gravity(v):
    return dot(f_vol, v) * dx


# Update functions

# Update acceleration
def update_a(u, u_old, v_old, a_old, ufl=True):
    if ufl:
        dt_ = dt
        beta_ = beta
    else:
        dt_ = float(dt)
        beta_ = float(beta)

    return ((u - u_old - dt_ * v_old) / beta / dt_ ** 2
            - (1 - 2 * beta_) / 2 / beta_ * a_old)


# Update velocity
def update_v(a, u_old, v_old, a_old, ufl=True):
    if ufl:
        dt_ = dt
        gamma_ = gamma
    else:
        dt_ = float(dt)
        gamma_ = float(gamma)

    return v_old + dt_ * ((1 - gamma_) * a_old + gamma_ * a)


def update_fields(u, u_old, v_old, a_old):
    """Update all fields at the end of a timestep."""

    u_vec, u0_vec = u.vector(), u_old.vector()
    v0_vec, a0_vec = v_old.vector(), a_old.vector()

    # call update functions
    a_vec = update_a(u_vec, u0_vec, v0_vec, a0_vec, ufl=False)
    v_vec = update_v(a_vec, u0_vec, v0_vec, a0_vec, ufl=False)

    # assign u->u_old
    v_old.vector()[:], a_old.vector()[:] = v_vec, a_vec
    u_old.vector()[:] = u.vector()


def avg(x_old, x_new, alpha):
    return alpha * x_old + (1 - alpha) * x_new


# residual
a_np1 = update_a(du, u_n, v_n, a_n, ufl=True)
v_np1 = update_v(a_np1, u_n, v_n, a_n, ufl=True)

res = m(avg(a_n, a_np1, alpha_m), v) + k(avg(u_n, du, alpha_f), v) # - W_gravity(v)

a_form = lhs(res)
L_form = rhs(res)

# parameters for Time-Stepping
t = 0.0
n = 0
E_ext = 0

write_interval = 0.05  # Salva ogni x secondi fisici
last_write_time = 0.0 # Tiene traccia dell'ultima scrittura

t_start_total = time.time()

displacement_out = File("output/u_NASA_g.pvd")

u_n.rename("Displacement", "")
u_np1.rename("Displacement", "")
displacement_out << (u_n, t)

print("\n[SOLID] Inizio ciclo di accoppiamento FSI...")

while precice.is_coupling_ongoing():

    if precice.requires_writing_checkpoint():  # write checkpoint
        # print(f"--- [SOLID] Salvataggio checkpoint al tempo t={t:.4f} ---")
        precice.store_checkpoint((u_n, v_n, a_n), t, n)
    
    # print(f"--- [SOLID] Attesa dati dal fluido (Time Window: {t:.4f}) ---")

    precice_dt = precice.get_max_time_step_size()
    dt = Constant(np.min([precice_dt, fenics_dt]))

    t_wait_start = time.time()
    read_data = precice.read_data((1 - float(alpha_f)) * dt)
    t_wait = time.time() - t_wait_start

    # read data from preCICE and get a new coupling expression
    # sample force F at $F(t_{n+1-\alpha_f})$ (see generalized alpha paper)
    read_data = precice.read_data((1 - float(alpha_f)) * dt)

    # print("--- [SOLID] Dati ricevuti. È il mio turno: risoluzione sistema strutturale ---")

    # Update the point sources on the coupling boundary with the new read data
    Forces_x, Forces_y, Forces_z = precice.get_point_sources(read_data)

    # print("      > Assemblo il sistema lineare (Matrice A e vettore b)...")
    A, b = assemble_system(a_form, L_form, bc)

    b_forces = b.copy()  # b is the same for every iteration, only forces change

    for ps in Forces_x:
        ps.apply(b_forces)
    for ps in Forces_y:
        ps.apply(b_forces)
    for ps in Forces_z:
        ps.apply(b_forces)

    assert (b is not b_forces)

    t_solve_start = time.time()

    # print(f"      > Risoluzione sistema lineare (Gradi di libertà: {V.dim()})...")
    solve(A, u_np1.vector(), b_forces)
    # print("      > Sistema risolto con successo.")

    # solver = PETScKrylovSolver("cg", "hypre_amg")
    # solver.solve(A, u_np1.vector(), b_forces)

    # solver = PETScKrylovSolver("gmres", "ilu")
    # solver.solve(A, u_np1.vector(), b_forces)
    
    t_solve = time.time() - t_solve_start

    t_precice_start = time.time()

    # Write new displacements to preCICE
    precice.write_data(u_np1)

    # print("--- [SOLID] Invio spostamenti a preCICE e calcolo convergenza ---")
    # Call to advance coupling, also returns the optimum time step value
    precice.advance(float(dt))

    t_precice = time.time() - t_precice_start

    # Either revert to old step if timestep has not converged or move to next timestep
    if precice.requires_reading_checkpoint():  # roll back to checkpoint
        # print("--- [SOLID] Convergenza NON raggiunta: ritorno al checkpoint (sotto-iterazione) ---")
        uva_cp, t_cp, n_cp = precice.retrieve_checkpoint()
        u_cp, v_cp, a_cp = uva_cp
        u_n.assign(u_cp)
        v_n.assign(v_cp)
        a_n.assign(a_cp)
        t = t_cp
        n = n_cp
    else:
        update_fields(u_np1, u_n, v_n, a_n)
        u_n.assign(u_np1)
        t += float(dt)
        n += 1

    if precice.is_time_window_complete():
        print(
        f"[SOLID] [TIMING] "
        f"wait={t_wait:.3f}s | "
        f"solve={t_solve:.3f}s | "
        f"precice={t_precice:.3f}s | "
        f"t={t:.3f}"
    )
        if t >= last_write_time + write_interval:
            displacement_out << (u_n, t)
            last_write_time = t
            print(f">>> [SOLID] Output salvato su file al tempo {t:.3f} <<<")

# Plot tip displacement evolution
displacement_out << (u_n, t)

t_end_total = time.time()
print(f"[SOLID] ClockTime totale = {t_end_total - t_start_total:.1f} s")

print("\n[SOLID] Simulazione completata. Finalizzazione...")
precice.finalize()
