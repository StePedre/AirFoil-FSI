# 3D Fluid-Structure Interaction of NASA SC(2)-0410 Airfoil

**Team:** Mattia Gotti, Michele Milani, Stefano Pedretti, Federico Pinto

## 📌 Overview
This project was developed for the High Performance Simulation Lab for Mechanical Engineering course. It presents a full three-dimensional (3D) Fluid-Structure Interaction (FSI) simulation of a deformable wing based on the NASA SC(2)-0410 airfoil. The simulation employs a high-fidelity partitioned approach to manage the complex aerodynamic and structural coupling implicitly.

## 🛠️ Technologies
* **Fluid Solver:** OpenFOAM
* **Solid Solver:** FEniCS
* **Coupling Library:** preCICE
* **Core Concepts:** 3D Fluid-Structure Interaction (FSI), implicit coupling, aerodynamic scalability analysis.

## 🚀 Key Features
* **Multi-Solver Integration:** Seamlessly couples OpenFOAM for the fluid domain and FEniCS for the structural domain using the preCICE open-source library.
* **High-Fidelity 3D Simulation:** Transitions from traditional 2D studies to a full 3D configuration to accurately capture complex, real-world aerodynamic interactions.
* **Scalability & Performance:** Includes a comprehensive scalability analysis across multiple cores (up to 16) to evaluate parallel efficiency and speedup under varying Angles of Attack (AoA) and coarse/refined meshes.

---
