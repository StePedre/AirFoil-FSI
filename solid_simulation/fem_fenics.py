from fenics import *
import numpy as np

# --- 1. PARAMETRI FISICI ---
E = 71.7e9    # Modulo di Young (E)
nu = 0.33            # Coefficiente di Poisson (nu)
rho = 2810.0         # Densità (rho)
g = 9.8             # Gravità (m/s^2)
TOL = 1E-10         # Tolleranza geometrica
dim = 3

# Parametri Materiale per il Report
MATERIALE = "Alluminio 7075-T6"

# Costanti di Lamé
mu = Constant(E / (2.0 * (1.0 + nu)))
lambda_ = Constant(E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu)))

# -------------------------------------------------------------------------
# FUNZIONE MAIN (Per replicare la struttura di output di test_fem.py)
# -------------------------------------------------------------------------
def main():
    
    # --- 2. SETUP INIZIALE E MESH ---
    print("\n" + "="*60)
    print(f"{'ANALISI ELASTICA LINEARE (FEM - FEniCS)':^60}")
    print("="*60)
    print(f"{'Materiale':<30} : {MATERIALE}")
    print(f"{'Modulo di Young (E)':<30} : {E/1e9:.1f} GPa")
    print(f"{'Coefficiente di Poisson (nu)':<30} : {nu:.2f}")
    print(f"{'Densità (rho)':<30} : {rho:.1f} kg/m³")
    print(f"{'Gravità (g)':<30} : {g:.2f} m/s²")
    
    print(" > Caricamento mesh...")

    # 1. Importa la Mesh (dominio 3D)
    mesh = Mesh()
    with XDMFFile("../mesh/mesh_domain.xdmf") as infile:
        infile.read(mesh)

    # Trova le coordinate estreme per definire geometricamente i contorni
    z_coordinates = mesh.coordinates()[:, 2] # Colonna Z
    Z_MIN = np.min(z_coordinates)
    Z_MAX = np.max(z_coordinates)

    print(f"   [Geo Check] Z_MIN trovato: {Z_MIN:.4f}")
    print(f"   [Geo Check] Z_MAX trovato: {Z_MAX:.4f}")

    # Spazio funzionale
    V = VectorFunctionSpace(mesh, 'P', 2)
    
    # --- 3. CONDIZIONI AL CONTORNO (BC) E FORMULAZIONE ---
    print(" > Assemblaggio sistema variazionale...")

    def clamped_boundary(x, on_boundary):
        # Incastro su Z_MIN
        return on_boundary and abs(x[2] - Z_MIN) < TOL 

    # Incastro (Dirichlet BC: u=0)
    fixed_boundary = AutoSubDomain(clamped_boundary)
    bc = DirichletBC(V, Constant((0.0, 0.0, 0.0)), fixed_boundary)

    u = TrialFunction(V)
    v = TestFunction(V)

    # Funzioni di Tensione e Deformazione
    def epsilon(u):
        return sym(grad(u))

    def sigma(u):
        return lambda_ * div(u) * Identity(dim) + 2 * mu * epsilon(u)

    # Carico di volume (Peso Proprio: g in direzione Y negativa)
    f_vol = Constant((0.0, -rho * g, 0.0))

    # Equazione di Navier per Elastostatica
    a = inner(sigma(u), epsilon(v)) * dx  # Lato Sinistro (Rigidità)
    L = dot(f_vol, v) * dx               # Lato Destro (Carico)

    # --- 4. RISOLUZIONE ---
    print(" > Risoluzione sistema lineare...")

    u_sol = Function(V)
    # L'output "Solving linear variational problem" viene soppresso dalla stampa
    solve(a == L, u_sol, bc) 
    u_sol.rename("Displacement", "")

    # --- 5. POST-PROCESSING E REPORT ---
    print(" > Calcolo stress Von Mises e altri scalari...")

    # Spazio discontinuo di grado 0 (DG0) per gli sforzi
    W = FunctionSpace(mesh, 'DG', 0) 
    sigma_u = sigma(u_sol)

    # TENSore DEVIATORICO e SFORZO DI VON MISES
    s = dev(sigma_u) 
    von_Mises_expr = sqrt(1.5 * inner(s, s))
    # NOTA: project() può generare output JIT se chiamato per la prima volta
    von_mises_sol = project(von_Mises_expr, W) 
    von_mises_sol.rename("VonMises_Stress", "")

    # SIGMA_ZZ 
    sigma_zz_expr = sigma_u[2, 2]
    sigma_zz_sol = project(sigma_zz_expr, W)
    sigma_zz_sol.rename("Sigma_zz", "")

    # SFORZO IDROSTATICO
    hydro_stress_expr = (1/3.0) * tr(sigma_u)
    hydro_stress_sol = project(hydro_stress_expr, W)
    hydro_stress_sol.rename("Hydrostatic_Stress", "")

    # Calcolo risultati massimi per il report
    u_array = u_sol.vector().get_local()
    u_magnitude = np.sqrt(u_array[::3]**2 + u_array[1::3]**2 + u_array[2::3]**2)
    max_displacement_m = np.max(u_magnitude)

    max_von_mises = von_mises_sol.vector().max()
    max_abs_sigma_zz = np.max(np.abs(sigma_zz_sol.vector().get_local()))
    max_abs_hydro_stress = np.max(np.abs(hydro_stress_sol.vector().get_local()))

    # --- REPORT FINALE UNIFORME ---
    print("\n" + "="*60)
    print(f"{'RISULTATI DEL CALCOLO':^60}")
    print("="*60)
    print(f"{'Spostamento Max (Tip Defl.)':<30} : {max_displacement_m * 1000.0:.4f} mm")
    print(f"{'Stress Von Mises Max':<30} : {max_von_mises / 1e6:.2f} MPa")
    print(f"{'Stress Sigma_zz Max (Modulo)':<30} : {max_abs_sigma_zz / 1e6:.2f} MPa")
    print(f"{'Stress Idrostatico Max (Modulo)':<30} : {max_abs_hydro_stress / 1e6:.2f} MPa")
    print("-" * 60)
    
    # --- 6. SALVATAGGIO OUTPUT ---
    print("\n > Salvataggio file XDMF unificato...")

    # Spostamento Vettoriale
    file_u = XDMFFile("results_fenics/displacement_fenics.xdmf")
    file_u.write(u_sol)
    file_u.close()

    # Sforzi Scalari (Usiamo PVD per raggruppare i campi scalari in FEniCS 2019)
    file_stress_pvd = File("results_fenics/stress_scalars_fenics.pvd")
    file_stress_pvd << von_mises_sol
    file_stress_pvd << sigma_zz_sol
    file_stress_pvd << hydro_stress_sol


    print(" > File di risultato salvati nella cartella results:")
    print("     • displacement_fenics.xdmf (Spostamento Vettoriale)")
    print("     • stress_scalars_fenics.pvd (Sforzi Scalari: Von Mises, Sigma_zz, Idrostatico)")
    print("\n")

if __name__ == "__main__":
    main()