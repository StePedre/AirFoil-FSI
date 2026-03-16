import os
import shutil
import argparse
import sys

def setup_simulation(angle, mesh_type, cpu_count):
    # ==========================================
    # 1. CONFIGURATION OF RELATIVE PATHS
    # ==========================================

    # Current folder
    current_case_path = os.getcwd()

    # Parent folder
    root_run_path = os.path.abspath(os.path.join(current_case_path, ".."))

    # Library paths
    meshes_library_path = os.path.join(root_run_path, "meshes")
    stls_library_path = os.path.join(root_run_path, "stlFiles")
    cpus_library_path = os.path.join(root_run_path, "cpuFiles")

    # ==========================================
    # 2. STL SELECTION LOGIC
    # ==========================================

    stl_map = {
        0:  "NASAsc2-0410_singleLine.stl",
        5:  "NASAsc2-0410_singleLine_5deg.stl",
        10: "NASAsc2-0410_singleLine_10deg.stl"
    }

    source_stl_filename = stl_map[angle]
    src_stl_path = os.path.join(stls_library_path, source_stl_filename)
    dst_stl_path = os.path.join(current_case_path, "constant", "triSurface", "NASAsc2-0410.stl")

    # ==========================================
    # 3. MESH DICT SELECTION LOGIC
    # ==========================================
    
    preset_folder_name = f"{angle}{mesh_type}"
    src_dict_path = os.path.join(meshes_library_path, preset_folder_name, "snappyHexMeshDict")
    dst_dict_path = os.path.join(current_case_path, "system", "snappyHexMeshDict")

    # ==========================================
    # 4. CPU FILE SELECTION LOGIC (Decompose & Jobs)
    # ==========================================

    # Specific folder for core count (e.g., ../cpuFiles/16)
    cpu_source_folder = os.path.join(cpus_library_path, str(cpu_count))

    # File A: decomposeParDict
    src_decomp_path = os.path.join(cpu_source_folder, "decomposeParDict")
    dst_decomp_path = os.path.join(current_case_path, "system", "decomposeParDict")

    # Files B: Job Scripts (List of scripts to copy)
    job_filenames = [
        "jobPimpleFoamJob_full.sh",
        "jobPimpleFoamJob_sim.sh",
        "jobPimpleFoamJob_mesh.sh"
    ]

    # ==========================================
    # 5. SAFETY CHECKS
    # ==========================================

    print(f"--- Case Setup: Angle {angle}° | Mesh {mesh_type} | CPU {cpu_count} ---")

    # Check STL existence
    if not os.path.exists(src_stl_path):
        print(f"\n[ERROR] STL file not found: {src_stl_path}")
        print(f"Check folder '../stlFiles/'")
        sys.exit(1)

    # Check Mesh Dict existence
    if not os.path.exists(src_dict_path):
        print(f"\n[ERROR] snappyHexMeshDict not found: {src_dict_path}")
        print(f"Check folder '../meshes/{preset_folder_name}/'")
        sys.exit(1)

    # Check CPU folder existence
    if not os.path.exists(cpu_source_folder):
        print(f"\n[ERROR] CPU folder not found: {cpu_source_folder}")
        print(f"Check existence of '../cpuFiles/{cpu_count}/'")
        sys.exit(1)

    # Check decomposeParDict existence
    if not os.path.exists(src_decomp_path):
        print(f"\n[ERROR] decomposeParDict missing in: {src_decomp_path}")
        sys.exit(1)

    # Check Job Scripts existence
    for job_file in job_filenames:
        src_job = os.path.join(cpu_source_folder, job_file)
        if not os.path.exists(src_job):
            print(f"\n[ERROR] Job script missing: {src_job}")
            sys.exit(1)

    # Create destination folders if missing
    os.makedirs(os.path.dirname(dst_stl_path), exist_ok=True)
    os.makedirs(os.path.dirname(dst_dict_path), exist_ok=True)

    # ==========================================
    # 6. COPY EXECUTION
    # ==========================================

    try:
        # Copy STL
        shutil.copy2(src_stl_path, dst_stl_path)
        print(f"✅ stl:  .../{source_stl_filename} -> constant/triSurface/")

        # Copy SnappyHexMeshDict
        shutil.copy2(src_dict_path, dst_dict_path)
        print(f"✅ mesh: .../{preset_folder_name}/snappyHexMeshDict -> system/")

        # Copy DecomposeParDict
        shutil.copy2(src_decomp_path, dst_decomp_path)
        print(f"✅ decomposeParDict:  .../cpuFiles/{cpu_count}/decomposeParDict -> system/")

        # Copy Job Scripts (Loop)
        for job_file in job_filenames:
            src_job = os.path.join(cpu_source_folder, job_file)
            dst_job = os.path.join(current_case_path, job_file)
            shutil.copy2(src_job, dst_job)
            print(f"✅ job:  .../cpuFiles/{cpu_count}/{job_file} -> ./{job_file}")

        print("\n--> Setup completed successfully.")

    except Exception as e:
        print(f"\n[CRITICAL] File copy error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Configure the case by copying STL, MeshDict, and CPU files.")

    # Flag Angle
    parser.add_argument(
        "-a", "--angle",
        type=int,
        choices=[0, 5, 10],
        required=True,
        help="Angle of attack (0, 5, 10)"
    )

    # Flag Mesh Type
    parser.add_argument(
        "-m", "--mesh",
        type=str,
        choices=["coarse", "refined"],
        required=True,
        help="Mesh type (coarse, refined)"
    )

    # Flag CPU
    parser.add_argument(
        "-c", "--cpu",
        type=int,
        choices=[1, 2, 4, 8, 16, 32],
        required=True,
        help="Number of CPUs (1, 2, 4, 8, 16, 32)"
    )

    args = parser.parse_args()

    setup_simulation(args.angle, args.mesh, args.cpu)
