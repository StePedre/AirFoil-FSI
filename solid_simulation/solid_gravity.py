# Import required libs
from fenics import Constant, Function, AutoSubDomain, VectorFunctionSpace, TrialFunction, TestFunction, DirichletBC, Identity, inner, dx, sym, grad, div, lhs, rhs, File, solve, assemble_system, Mesh, XDMFFile, dot
import numpy as np

# 1. Importa la Mesh (dominio 3D)
mesh = Mesh()
with XDMFFile("mesh/mesh_domain_coarse.xdmf") as infile:
    infile.read(mesh)


# --- PARAMETRI FISICI e GEOMETRICI ---
dim = 3
E = 71.7e9          # Modulo di Young (E)
nu = 0.33           # Coefficiente di Poisson (nu)
rho = 2810.0        # Densità (rho)

g = 9.81            # Gravità (m/s^2)
f_vol = Constant((0.0, -rho * g, 0.0))

mu = Constant(E / (2.0 * (1.0 + nu)))
lambda_ = Constant(E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu)))

# --- SPAZIO FUNZIONALE E FUNZIONI ---
V = VectorFunctionSpace(mesh, 'P', 1)

# Trial and Test Functions
du = TrialFunction(V)
v = TestFunction(V)

u_np1 = Function(V) # Soluzione al passo n+1
u_n = Function(V)   # Spostamento al passo n
v_n = Function(V)   # Velocità al passo n
a_n = Function(V)   # Accelerazione al passo n

# Inizializzazione (a zero)
u_n.vector()[:] = 0.0
v_n.vector()[:] = 0.0
a_n.vector()[:] = 0.0

# --- CONDIZIONI AL CONTORNO ---
TOL = 1E-10
z_coordinates = mesh.coordinates()[:, 2] 
Z_MIN = np.min(z_coordinates)

def clamped_boundary(x, on_boundary):
    # Incastro sulla superficie Z=Z_MIN
    return on_boundary and abs(x[2] - Z_MIN) < TOL 
    
fixed_boundary = AutoSubDomain(clamped_boundary)
bc = DirichletBC(V, Constant((0.0, 0.0, 0.0)), fixed_boundary)


# --- PARAMETRI DI INTEGRAZIONE TEMPORALE ---

alpha_m = Constant(0.2) 
alpha_f = Constant(0.4) 

gamma = Constant(0.5 + alpha_f - alpha_m)
beta = Constant((gamma + 0.5)**2 / 4.)

# Passo temporale FENICS (Fisso, non accoppiato)
DT = 0.01
T_END = 5.0
dt = Constant(DT)

# --- DEFINIZIONE FORME (LE STESSE DI solid.py) ---

def epsilon(u):
    return 0.5 * (grad(u) + grad(u).T)

def sigma(u):
    return lambda_ * div(u) * Identity(dim) + 2 * mu * epsilon(u)

def m(u, v):
    return rho * inner(u, v) * dx

def k(u, v):
    return inner(sigma(u), sym(grad(v))) * dx

def W_gravity(v):
    return dot(f_vol, v) * dx

def avg(x_old, x_new, alpha):
    return alpha * x_old + (1 - alpha) * x_new

# --- FUNZIONI DI AGGIORNAMENTO DEI CAMPI (Le stesse di solid.py) ---

def update_a(u, u_old, v_old, a_old, ufl=True):
    # ... (Il corpo della funzione è identico a solid.py)
    if ufl:
        dt_ = dt
        beta_ = beta
    else:
        dt_ = float(dt)
        beta_ = float(beta)
    return ((u - u_old - dt_ * v_old) / beta / dt_ ** 2 - (1 - 2 * beta_) / 2 / beta_ * a_old)

def update_v(a, u_old, v_old, a_old, ufl=True):
    # ... (Il corpo della funzione è identico a solid.py)
    if ufl:
        dt_ = dt
        gamma_ = gamma
    else:
        dt_ = float(dt)
        gamma_ = float(gamma)
    return v_old + dt_ * ((1 - gamma_) * a_old + gamma_ * a)

def update_fields(u, u_old, v_old, a_old):
    """Aggiorna tutti i campi alla fine di un passo temporale."""
    u_vec, u0_vec = u.vector(), u_old.vector()
    v0_vec, a0_vec = v_old.vector(), a_old.vector()
    
    # Calcola e assegna i nuovi valori
    a_vec = update_a(u_vec, u0_vec, v0_vec, a0_vec, ufl=False)
    v_vec = update_v(a_vec, u0_vec, v0_vec, a0_vec, ufl=False)

    v_old.vector()[:], a_old.vector()[:] = v_vec, a_vec
    u_old.vector()[:] = u.vector()

# --- FORMULAZIONE DEL SISTEMA LINEARE ---

# Termini n+1 (TrialFunction)
a_np1 = update_a(du, u_n, v_n, a_n, ufl=True)

# Il Residuale include la forza peso W_gravity, ma NON le forze FSI
res = m(avg(a_n, a_np1, alpha_m), v) + k(avg(u_n, du, alpha_f), v) - W_gravity(v)

a_form = lhs(res)
L_form = rhs(res)

# Assembla la matrice A e il vettore b solo una volta se sono statici (qui b non lo è!)
# Poiché L_form contiene termini u_n, v_n, a_n, A e b devono essere riassemblati
# in ogni passo (o ricompilati UFL/JIT se FEniCS lo permette).

# --- CICLO TEMPORALE AUTONOMO ---
t = 0.0
n = 0
displacement_out = File("output_gravity/u_gravity.pvd")

u_n.rename("Displacement", "")
displacement_out << (u_n, t) # Salva la condizione iniziale

# Report Iniziale
print("\n" + "="*50)
print(f"{'SIMULAZIONE DINAMICA - SOLO GRAVITÀ':^50}")
print("="*50)
print(f"{'Tempo Finale (T_END)':<30} : {T_END:.3f} s")
print(f"{'Passo Temporale (DT)':<30} : {DT:.5f} s")
print(f"{'Numero di Nodi Mesh':<30} : {mesh.num_vertices()}")
print(f"{'Direzione Gravità':<30} : Y negativa")
print("-" * 50)
print(f"{'PASSAGGIO':<10}{'TEMPO (s)':<15}{'MAX DISPL (mm)':<15}")
print("-" * 50)

while t < T_END:
    
    # Risoluzione del sistema:
    # A e b vengono assemblati, b contiene implicitamente il carico di gravità
    A, b = assemble_system(a_form, L_form, bc)

    # In questo script, NON c'è bisogno di copiare b e applicare Forces_x/y/z.
    # Il vettore b è il carico totale (che qui è solo la Gravità).
    solve(A, u_np1.vector(), b)

    # Aggiornamento dei campi per il passo successivo:
    update_fields(u_np1, u_n, v_n, a_n)
    u_n.assign(u_np1)
    
    # Avanzamento del tempo:
    t += float(dt)
    n += 1

    if n % 10 == 0:
        # Calcolo del modulo del vettore spostamento u_n
        u_array = u_n.vector().get_local()
        # V è un VectorFunctionSpace di dimensione 3 (u_x, u_y, u_z)
        u_magnitude = np.sqrt(u_array[::3]**2 + u_array[1::3]**2 + u_array[2::3]**2)
        max_displacement_m = np.max(u_magnitude)
        max_displacement_mm = max_displacement_m * 1000.0

        print(f"{n:<10}{t:<15.3f}{max_displacement_mm:<15.4f}")
        displacement_out << (u_n, t)

# Salvataggio finale e Report
displacement_out << (u_n, t)
print("-" * 50)
print(f"Simulazione completata. Tempo finale raggiunto: {t:.3f} s")
print(f"Risultati salvati in output_gravity/u_gravity.pvd")