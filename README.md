# CNDO/2 Self-Consistent Field Program

## University of California, Berkeley  
**CHEM279: Numerical Algorithms Applied to Computational Quantum Chemistry**  
**SCF Program - CNDO/2 Method**  
**Author: Brandon Robello**

---

## Objective:
This repository contains a CNDO/2 (Complete Neglect of Differential Overlap) **self-consistent field (SCF)** program to evaluate molecular orbitals and total energy of molecules with open or closed shells. The project advances previous work from extended Huckel theory, leveraging the semi-empirical Hartree-Fock method, which provides both computational efficiency and conceptual simplicity.

This program constructs basis functions from **Gaussian-type orbitals (GTOs)**, assembles matrices (such as overlap, Fock, and Hamiltonian matrices), and iteratively solves the SCF equations for molecular orbitals and energy.

The program applies the CNDO/2 method to:

1. **H2** (closed-shell molecule)
2. **HO•** (open-shell molecule with 4 alpha electrons and 3 beta electrons)
2. **HF** (closed-shell molecule)
2. **N2** (closed-shell molecule)
2. **N** (open-shell molecule with 4 alpha electrons and 1 beta electrons)
---

## Project Goals:
1. **Build and Assemble the CNDO/2 Fock Matrix**:
   - Implement functions to calculate **diagonal and off-diagonal Fock matrix elements** for alpha and beta spin states.
   - Use semi-empirical parameters for H, C, N, O, and F atoms provided in the problem statement, ensuring consistency in units (atomic units to electron volts).
   
2. **Self-Consistent Field (SCF) Procedure**:
   - Start with an initial guess for electron density matrices (Initial Guess of 0 for this program).
   - Build Fock matrices and solve the generalized eigenvalue problem to obtain molecular orbital coefficients and energies.
   - Iterate until convergence, **updating the density matrices with each iteration**.

3. **Total Energy Calculation**:
   - Calculate the **total electronic energy** from the converged SCF results.

4. **Exploratory Chemistry**:
   - Extend the application of CNDO/2 to molecules to N2 to calculate bond energy.

---

## File Structure:

```
|--- CMakeLists.txt               # CMake configuration
|--- Include/
|   |--- Atom.h                   # Atom class declaration
|   |--- Molecule.h               # Molecule class declaration
|   |--- BasisFunction.h          # Basis function class declaration
|   |--- Reader.h                 # Reader abstract class
|   |--- AtomReader.h             # Reader for atom data
|--- Source/
|   |--- Atom.cpp                 # Atom class implementation
|   |--- Molecule_hartreeFock.cpp # Molecule class Hartree-Fock definitions
|   |--- Molecule_init.cpp        # Molecule class initialization definitions
|   |--- Molecule_overlap.cpp     # Molecule class Orbital Overlap definitions
|   |--- Molecule_utils.cpp     # Molecule class Utility definitions
|--- Tests/
|   |--- SCF_test.cpp           # Test executable for SCF CNDO/2 calculations
|--- Utils/
|   |--- Inputs/                  # Sample input files
|   |--- Outputs/                 # Output log files
```

---

## Running the Project

### Clone the Repository:
```bash
git clone https://github.com/brandonrobello/CHEM279_HW4.git
cd CHEM279_HW4
```

### Build the Project:
```bash
mkdir build
cd build
cmake ..
make
```

### Run the Test Case:
```bash
./SCF_test
```

### Output:
The output file for the test cases will be located in `Utils/Outputs/SCF_test_log.txt`.

---

## Key Modules and Calculations

### 1. Molecule and Basis Function Setup

The program reads molecular geometry in the format:
```
E X Y Z
```
Where `E` is the atomic symbol (e.g., H or O), and `X`, `Y`, `Z` are the Cartesian coordinates. The program constructs basis functions using **STO-3G** parameters (from Homework 4) for each atom's valence orbitals, initializing electron densities for alpha and beta spins as required for open-shell molecules.

### 2. Overlap Matrix Calculation

The **overlap matrix** $S_{\mu \nu}$ is calculated by summing over integrals of primitive Gaussian functions. For open-shell systems, this matrix serves as a key component of the SCF eigenvalue problem.

### 3. Fock Matrix Assembly

The **Fock matrix** $f_{\mu \nu}$ consists of diagonal and off-diagonal elements that depend on the density matrices, semi-empirical parameters, and overlap integrals:

- **Diagonal elements**: Contain terms for ionization energy, electron affinity, and electron repulsion ($\gamma$ terms).
- **Off-diagonal elements**: Depend on bonding parameters and overlap integrals.

### 4. SCF Iterative Solution

The **SCF algorithm** is implemented as follows:

1. **Guess Initial Density**: Set initial density matrices for alpha and beta electrons to zero.
2. **Build Fock Matrix**: Construct Fock matrices based on the current density matrices.
3. **Solve Eigenvalue Problem**: Obtain molecular orbital coefficients and energies.
4. **Update Density Matrices**: Recompute densities using the p lowest-energy alpha and q lowest-energy beta orbitals.
5. **Check for Convergence**: Iterate until changes in density matrices are below a specified threshold.

### 5. Total Energy Calculation

The **total electronic energy** is calculated from the converged SCF results as:
$$
E_{CNDO/2} = \frac{1}{2} \sum_{\mu \nu} p_{\mu \nu} (h_{\mu \nu} + f_{\mu \nu}) + E_{\text{nuclear}}
$$
where $E_{\text{nuclear}}$ is the nuclear repulsion energy calculated from atomic positions and charges.