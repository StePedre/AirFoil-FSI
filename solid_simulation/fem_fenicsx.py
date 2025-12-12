import numpy as np
import ufl
from dolfinx import fem, io, default_scalar_type, log
from dolfinx.fem.petsc import LinearProblem
from mpi4py import MPI
from petsc4py import PETSc

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

def main():
    # -------------------------------------------------------------------------
    # 1. CONFIGURAZIONE PROBLEMA
    # -------------------------------------------------------------------------
    E = 71.7e9          # 71.7 GPa
    nu = 0.33
    rho = 2810.0        
    g = 9.8   
    ID_VINCOLO = 8      # numero ID vincolo incastro
    MATERIALE = "Alluminio 7075-T6"

    mu = E / (2.0 * (1.0 + nu))
    lambda_ = (E * nu) / ((1.0 + nu) * (1.0 - 2.0 * nu))

    if rank == 0:
        print("\n" + "="*60)
        print(f"{'ANALISI ELASTICA LINEARE (FEM - DOLFINx)':^60}")
        print("="*60)
        print(f"{'Materiale':<30} : {MATERIALE}")
        print(f"{'Modulo di Young (E)':<30} : {E/1e9:.1f} GPa")
        print(f"{'Coefficiente di Poisson (nu)':<30} : {nu:.2f}")
        print(f"{'Densità (rho)':<30} : {rho:.1f} kg/m³")
        print(f"{'Gravità (g)':<30} : {g:.2f} m/s²")


    # -------------------------------------------------------------------------
    # 2. MESH
    # -------------------------------------------------------------------------
    if rank == 0: print(" > Caricamento mesh e conversione unità...")
    
    # Percorsi assunti per coerenza
    with io.XDMFFile(comm, "mesh_fenicsx/mesh_domain.xdmf", "r") as xdmf: 
        domain = xdmf.read_mesh(name="Grid")

    tdim = domain.topology.dim
    fdim = tdim - 1
    domain.topology.create_connectivity(tdim, fdim)
    domain.topology.create_connectivity(fdim, tdim)

    with io.XDMFFile(comm, "mesh_fenicsx/mesh_boundary.xdmf", "r") as xdmf: 
        boundary_tags = xdmf.read_meshtags(domain, name="Grid")

    ds = ufl.Measure("ds", domain=domain, subdomain_data=boundary_tags)
    dx = ufl.Measure("dx", domain=domain)

    # -------------------------------------------------------------------------
    # 3. DEFINIZIONE PROBLEMA
    # -------------------------------------------------------------------------
    if rank == 0: print(" > Assemblaggio sistema variazionale...")
    
    # Spazio funzionale P2
    V = fem.functionspace(domain, ("Lagrange", 2, (domain.geometry.dim,)))
    
    root_dofs = fem.locate_dofs_topological(V, fdim, boundary_tags.find(ID_VINCOLO))
    u_zero = np.array([0, 0, 0], dtype=default_scalar_type)
    bc_root = fem.dirichletbc(u_zero, root_dofs, V)

    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)

    def epsilon(u): return ufl.sym(ufl.grad(u))
    def sigma(u): return lambda_ * ufl.nabla_div(u) * ufl.Identity(len(u)) + 2 * mu * epsilon(u)

    f_vol = fem.Constant(domain, default_scalar_type((0, -g * rho, 0)))
    
    a = ufl.inner(sigma(u), epsilon(v)) * dx
    L = ufl.dot(f_vol, v) * dx

    # -------------------------------------------------------------------------
    # 4. SOLUZIONE
    # -------------------------------------------------------------------------
    if rank == 0: print(" > Risoluzione sistema lineare (MUMPS)...")
    
    problem = LinearProblem(a, L, bcs=[bc_root], 
                            petsc_options={"ksp_type": "preonly", "pc_type": "lu"},
                            petsc_options_prefix="elasticity_solve")
    u_sol = problem.solve()
    u_sol.name = "Displacement"

    # -------------------------------------------------------------------------
    # 5. POST-PROCESSING 
    # -------------------------------------------------------------------------
    W = fem.functionspace(domain, ("DG", 0))
    sigma_u = sigma(u_sol)
    
    # --- SFORZO DI VON MISES (Scalare) ---
    if rank == 0: print(" > Calcolo stress Von Mises e altri scalari...")
    deviatoric_stress = sigma_u - (1./3) * ufl.tr(sigma_u) * ufl.Identity(len(u_sol))
    von_Mises_expr = ufl.sqrt(1.5 * ufl.inner(deviatoric_stress, deviatoric_stress))
    
    stress_expr = fem.Expression(von_Mises_expr, W.element.interpolation_points)
    von_mises_sol = fem.Function(W)
    von_mises_sol.interpolate(stress_expr)
    von_mises_sol.name = "VonMises_Stress"
    
    # --- Componenti Scalari Specifiche ---
    # 1. Componente zz dello Sforzo (Sigma_zz)
    sigma_zz_expr = fem.Expression(sigma_u[2, 2], W.element.interpolation_points)
    sigma_zz_sol = fem.Function(W)
    sigma_zz_sol.interpolate(sigma_zz_expr)
    sigma_zz_sol.name = "Sigma_zz"

    # 2. Sforzo Idrostatico (Hydrostatic Stress)
    hydro_stress_expr = fem.Expression((1/3.) * ufl.tr(sigma_u), W.element.interpolation_points)
    hydro_stress_sol = fem.Function(W)
    hydro_stress_sol.interpolate(hydro_stress_expr)
    hydro_stress_sol.name = "Hydrostatic_Stress"

    # -------------------------------------------------------------------------
    # 6. REPORT
    # -------------------------------------------------------------------------
    u_vals = u_sol.x.array.reshape((-1, 3))
    local_max_u = np.max(np.linalg.norm(u_vals, axis=1))
    global_max_u = comm.allreduce(local_max_u, op=MPI.MAX)

    stress_vals = von_mises_sol.x.array
    local_max_s = np.max(stress_vals)
    global_max_s = comm.allreduce(local_max_s, op=MPI.MAX)
    
    sigma_zz_vals = sigma_zz_sol.x.array
    local_max_abs_s_zz = np.max(np.abs(sigma_zz_vals))
    global_max_abs_s_zz = comm.allreduce(local_max_abs_s_zz, op=MPI.MAX)
    
    hydro_stress_vals = hydro_stress_sol.x.array
    local_max_abs_hydro = np.max(np.abs(hydro_stress_vals))
    global_max_abs_hydro = comm.allreduce(local_max_abs_hydro, op=MPI.MAX)


    if rank == 0:
        # Report identico a solid_gravity.py
        print("\n" + "="*60)
        print(f"{'RISULTATI DEL CALCOLO':^60}")
        print("="*60)
        print(f"{'Spostamento Max (Tip Defl.)':<30} : {global_max_u * 1000:.4f} mm")
        print(f"{'Stress Von Mises Max':<30} : {global_max_s / 1e6:.2f} MPa")
        print(f"{'Stress Sigma_zz Max (Modulo)':<30} : {global_max_abs_s_zz / 1e6:.2f} MPa")
        print(f"{'Stress Idrostatico Max (Modulo)':<30} : {global_max_abs_hydro / 1e6:.2f} MPa")
        print("-" * 60)
        

    # -------------------------------------------------------------------------
    # 7. SALVATAGGIO UNIFICATO 
    # -------------------------------------------------------------------------
    if rank == 0: print("\n > Salvataggio file XDMF unificato...")

    # 1. Spostamento (u) - Campo Vettoriale
    with io.XDMFFile(comm, "results_fenicsx/displacement_fenicsx.xdmf", "w") as xdmf: 
        xdmf.write_mesh(domain)
        xdmf.write_function(u_sol) 
        
    # 2. Dati Scalari (Stress)
    with io.XDMFFile(comm, "results_fenicsx/stress_fenicsx.xdmf", "w") as xdmf: 
        xdmf.write_mesh(domain)
        xdmf.write_function(von_mises_sol) 
        xdmf.write_function(sigma_zz_sol)  
        xdmf.write_function(hydro_stress_sol) 
        
    

    if rank == 0:
        print(f" > File di risultato salvati nella cartella results:")
        print(f"     • displacement_fenicsx.xdmf (Spostamento Vettoriale)")
        print(f"     • stress_fenicsx.xdmf (Sforzi Scalari: Von Mises, Sigma_zz, Idrostatico)")
        print("\n")


if __name__ == "__main__":
    main()