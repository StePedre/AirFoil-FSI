## Scripts

The `pimpleFoam` directory contains four Python scripts designed to manage the simulation workflow:

* **`setup_case.py`**
    Sets up the simulation case by specifying the Angle of Attack (AoA), mesh refinement level, and CPU count.
    ```bash
    python3 setup_case.py -a [0, 5, 10] -m [coarse, refined] -c [1, 2, 4, 8, 16, 32]
    ```

* **`info_setup.py`**
    Displays the details of the configuration currently loaded.
    ```bash
    python3 info_setup.py
    ```

* **`print_avg_cl_cd.py`**
    Calculates and prints the average values of the lift ($C_l$) and drag ($C_d$) coefficients.
    ```bash
    python3 print_avg_cl_cd.py postProcessing/forceCoeff1/0/coefficient.dat
    ```

* **`print_yPlus.py`**
    Prints the average $y^+$ value on the walls and specifically on the `NASAsc2-0410` patch.
    ```bash
    python3 print_yPlus.py postProcessing/yPlus_monitoring/0/yPlus.dat
    ```