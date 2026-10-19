# Galactic_ISM_and_Magnetic_fields
Codes done in two separate projects (one discontinued) under Prof. Luke Chamandy 

## Computational Physics Research Portfolio

This repository contains recreated computational models and numerical simulations under the supervision of Prof. Luke Chamandy at NISER Bhubaneswar.

Both Projects focus on the working of Galactic magnetic fields and the Galactic Dynamo. The focus of this work has been translating theoretical frameworks into functional algorithms and replicating peer-reviewed baselines to establish numerical accuracy.

---

## 🔬 Project 1: [Modelling Turbulence in the Galactic interstellar medium]
Status: Completed (Code Validation Framework)

This project focused on translating an established physical framework into a functional numerical solver by replicating an existing model from a student paper. 

*   Objective: To model how supernovae and superbubbles inject energy into the galactic interstellar medium and thus drive turbulence.
*   Methodology & Implementations:-
 - Implemented a spatially varying gas-density model incorporating atomic hydrogen (HI) and molecular hydrogen (H₂).
 - Modelled the spatial distribution of supernova rates and calculated the characteristic sizes, expansion velocities, and energy contributions of supernova   remnants and superbubbles.
 - Estimated the turbulent root-mean-square (RMS) velocity, integral length scale, eddy turnover time, and outflow velocity across the galactic disc and halo.
 - Investigated how the spatial distribution of energy injection and the relative contributions of isolated and clustered supernovae affect turbulence.
*   Key Outcome: Successfully recreated the reference model's key plots validating the codebase as a reliable benchmark for this physical system.

### 📖 Reference Literature -
 - Internal research report by Daniel Turnbull from the research group of Prof. Anvar Shukurov(Newcastle university).
 - Internal research report by Ryan Smith from the research group of Prof. Anvar Shukurov(Newcastle university).

 This code is a replication of the computational models developed in above papers provided via academic collaboration through Dr. Luke Chamandy. It serves as a benchmark suite to validate numerical accuracy.

---

## 🛰️ Project 2: [Modelling the Magnetic field of Galaxies]
Status: Ongoing (Active Development)

To investigate the generation and evolution of large-scale magnetic fields in the Milky Way through galactic dynamo theory, focusing on the interplay between magnetic-field components and turbulent processes.

*   Current Progress: Successfully translated and replicated the foundational layers of a larger theoretical model.
*   Methodology and Implementation:-
  - Study the mean-field dynamo framework, which describes how turbulent motions and differential rotation contribute to the generation and maintenance of large-scale galactic magnetic fields.
  - Model the coupled evolution of the radial magnetic field (B_r) and azimuthal magnetic field (B_phi), incorporating magnetic diffusion and dynamo cpling.
  - Investigate nonlinear magnetic-helicity feedback through the magnetic contribution to the dynamo alpha effect (alpha_m), which can suppress further magnetic-field growth.
  - Develop a computational approach involving spatial discretisation and Runge–Kutta time integration to study magnetic-field evolution and nonlinear saturation

### 📖 Reference Literature -
 - Chamandy, Luke, Kandaswamy Subramanian, and Anvar Shukurov. "Galactic spiral patterns and dynamo action–I. A new twist on magnetic arms." Monthly Notices of the Royal Astronomical Society 428.4 (2013): 3569-3589.

 This code is a replication of the computational model A developed in above paper.


##  Languages & Libraries Used
- Python 3
- math (standard library)
- NumPy,Scipy (for numerical methods, curvefitting etc)
- Matplotlib (for plotting and data visualization)
