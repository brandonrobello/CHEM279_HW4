// Chem279: Numerical Algorithms Applied to Computational Quantum Chemistry  
// Creator: Brandon Robello
// Date Adapted: 11/8/24
//
// Molecule.h contains the necessary include files and 
//        C++ Class declarations for Molecule class


#include "Atom.h"
#include "BasisFunction.h"
#include <armadillo>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
#include <set>
#include <map>

#ifndef MOLECULE_H
#define MOLECULE_H



// Define a struct to hold the semi-empirical parameters for each element
struct CNDO2Parameters {
    double IA_s;  ///< Semi-empirical parameter for s-orbital (1/2 (Is + As)).
    double IA_p;  ///< Semi-empirical parameter for p-orbital (1/2 (Ip + Ap)).
    double beta;  ///< Semi-empirical parameter for bonding interaction (-β).
};


/**
 * @class Molecule
 * 
 * @brief Represents a molecular system in quantum chemistry with its atoms, basis functions, and quantum matrices.
 * 
 * The Molecule class stores information about the atoms, basis functions, and quantum matrices such as
 * the overlap matrix, Hamiltonian matrix, and molecular orbital coefficients. It supports computing overlaps,
 * building the Hamiltonian, and calculating the molecule's total energy.
 * 
 * @details
 * - **atoms_**: Vector storing the atoms that make up the molecule.
 * - **basisFunctions_**: Vector of basis functions used to describe molecular orbitals.
 * - **S_overlap_matrix_**: Overlap matrix for the basis functions.
 * - **H_matrix_**: Hamiltonian matrix for the molecule.
 * - **C_matrix_**: Molecular orbital coefficient matrix.
 * 
 * @see Atom, BasisFunction
 */
class Molecule {
public:
    /**
     * @brief Constructs a Molecule object from a vector of atoms.
     * 
     * This constructor takes a vector of Atom objects and builds the molecule by initializing its atoms
     * and computing associated properties, such as basis functions.
     * 
     * @param Atoms A vector of Atom objects representing the atoms of the molecule.
     */
    Molecule(const std::vector<Atom>& Atoms);

    /**
     * @brief Overload for the << operator to print the Molecule details.
     * 
     * This function prints the relevant information of the molecule in a human-readable format,
     * including atoms and their properties.
     * 
     * @param os Output stream to print the molecule information.
     * @param molecule The Molecule object to be printed.
     * @return std::ostream& Reference to the output stream with molecule details.
     */
    friend std::ostream& operator<<(std::ostream& os, const Molecule& molecule);

    /**
     * @brief Adds a new atom to the molecule.
     * 
     * Adds an Atom object to the molecule and updates the basis functions.
     * 
     * @param atom The Atom object to be added.
     */
    void addAtom(const Atom& atom);

    /**
     * @brief Returns the vector of atoms in the molecule.
     * @return std::vector<Atom> A vector of atoms.
     */
    std::vector<Atom> getAtoms() const;

    /**
     * @brief Prints the basis functions of the molecule.
     */
    void printBasisFunctions() const;

    /**
     * @brief Prints the overlap matrix of the molecule.
     */
    void printS() const;

    /**
     * @brief Prints the electron repulsion integral (gamma) matrix of the molecule.
     */
    void printGamma() const;
    
    /**
     * @brief Prints the Hamiltonian (H matrix) of the molecule.
     */ 
    void printH() const;

    /**
     * @brief Calculates the contracted overlap between two basis functions.
     * 
     * This function computes the overlap between two contracted Gaussian-type orbitals (basis functions).
     * 
     * @param mu Index of the first basis function.
     * @param nu Index of the second basis function.
     * @return double The computed overlap.
     */
    double calculate_contracted_overlap(int mu, int nu) const;

    /**
     * @brief Calculates the primitive overlap between two Gaussian primitives.
     * 
     * Computes the overlap integral between two uncontracted Gaussian primitives based on their exponents
     * and positions.
     * 
     * @param alpha_k Exponent of the first Gaussian.
     * @param alpha_l Exponent of the second Gaussian.
     * @param R_A Center of the first Gaussian.
     * @param R_B Center of the second Gaussian.
     * @param lA Angular momentum of the first Gaussian.
     * @param lB Angular momentum of the second Gaussian.
     * @return double The computed primitive overlap integral.
     */
    double calculate_primitive_overlap(double alpha_k, double alpha_l, const arma::vec& R_A, const arma::vec& R_B, const arma::vec lA, const arma::vec lB) const;

    /**
     * @brief Constructs the overlap matrix for the molecule.
     */
    void calculate_overlap_matrix();

    /**
     * @brief Calculates the total electronic energy of the molecule.
     * 
     * Computes the total electronic energy of the molecule using the Hamiltonian and
     * molecular orbital coefficients.
     * 
     * @return double The total electronic energy of the molecule.
     */
    double calc_totalE();

    /**
     * @brief Builds the Fock matrix for the molecule with Self-Consistent Field (SCF) iterations.
     * 
     * @param verbose If true, prints additional debug information.
     */
    void build_Fmatrix_wSCF(bool verbose = false);

    /**
     * @brief Prints energies to a specified file.
     * 
     * @param file The filename to output energies to.
     */
    void printEnergies(const std::string& file) const;

private:
    /**
     * @brief Builds basis functions for a given atom.
     * 
     * Initializes the basis functions for an atom and adds them to the molecule's basis functions list.
     * 
     * @param atom The Atom object for which to build basis functions.
     */
    void buildBasisFunctions();
    
    /**
     * @brief Initializes molecule properties, including valence electrons and basis functions.
     */
    void Molecule_init();            

// Class members for storing atoms, basis functions, matrices, and properties
    std::vector<Atom> atoms_;                   ///< Vector containing the atoms of the molecule.
    std::vector<BasisFunction> basisFunctions_; ///< List of basis functions in the molecule.
    arma::mat S_overlap_matrix_;                ///< Overlap matrix for basis functions.
    arma::mat FA_matrix_;                       ///< Alpha Fock matrix.
    arma::mat FB_matrix_;                       ///< Beta Fock matrix.
    arma::mat pA_matrix_;                       ///< Alpha density matrix.
    arma::mat pB_matrix_;                       ///< Beta density matrix.
    arma::mat H_core_matrix_;                   ///< Core Hamiltonian matrix.
    arma::mat pTot_mu_nu_;                      ///< Total density matrix.
    arma::vec pTot_;                            ///< Total density vector for each atom.
    arma::mat VA_matrix_;                       ///< Alpha eigenvectors matrix.
    arma::mat VB_matrix_;                       ///< Beta eigenvectors matrix.
    arma::vec epsilon_A_;                       ///< Alpha eigenvalues vector.
    arma::vec epsilon_B_;                       ///< Beta eigenvalues vector.
    arma::mat cA_matrix_;                       ///< Alpha MO coefficients matrix.
    arma::mat cB_matrix_;                       ///< Beta MO coefficients matrix.
    arma::mat gamma_matrix_;                    ///< Gamma matrix for electron repulsion.
    int N_ = 0;                                 ///< Number of basis functions in the molecule.
    int val_electrons_ = 0;                     ///< Number of valence electrons in the molecule.
    int numCarbons_ = 0;                        ///< Number of carbon atoms.
    int numHydrogens_ = 0;                      ///< Number of hydrogen atoms.
    int alpha_e_, beta_e_ = 0;                  ///< Alpha and beta electron counts.
    double alpha_electron_E_, beta_electron_E_ = 0; ///< Energy associated with alpha and beta electrons.
    int charge_ = 0;                            ///< Total charge of the molecule.
    double nuclear_repulsion_E_ = 0;            ///< Nuclear repulsion energy.
    double total_energy_ = 0;                   ///< Total electronic energy.
    std::unordered_map<int, CNDO2Parameters> semiEmpiricalParams_ = { // Map with semi-empirical parameters for H, C, N, O, F
        {1, {7.176, 0.0, -9}},    // H: Ionization energy and beta (only s-orbital)
        {6, {14.051, 5.572, -21}}, // C: Ionization and p-orbital values
        {7, {19.316, 7.275, -25}}, // N: Ionization and p-orbital values
        {8, {25.390, 9.111, -31}}, // O: Ionization and p-orbital values
        {9, {32.272, 11.080, -39}} // F: Ionization and p-orbital values
    };

    // Helper Functions for build and calculate functions
    /**
     * @brief Calculates the exponential prefactor for Gaussian overlap in one dimension.
     */
    double exponential_prefactor(double alpha, double beta, double X_A, double X_B) const;
    /**
     * @brief Computes the factorial of a number (n!).
     */
    int factorial(int n) const;
    /**
     * @brief Computes the double factorial of a number (n!!).
     */
    int double_factorial(int n) const;
    /**
     * @brief Calculates the binomial coefficient \f$ \binom{l}{i} \f$.
     */
    double calculate_binomial(int l, int i) const;
    /**
     * @brief Calculates the center of the product Gaussian between two primitives.
     */
    double calculate_RP(double alpha, double beta, double X_A, double X_B) const;
    /**
     * @brief Computes the power term \f$ (R_p - X_A)^{l-i} \f$ used in overlap integrals.
     */
    double calculate_power_term(double Rp, double center, int l, int i) const;
    /**
     * @brief Computes the term used in the one-dimensional overlap matrix equation.
     */
    double calculate_term(double alpha, double beta, double X_A, double X_B, int lA, int lB, int i, int j) const;
    /**
     * @brief Calculates the one-dimensional overlap integral between two Gaussian primitives.
     */
    double Overlap_onedim(double alpha, double beta, double X_A, double X_B, int lA, int lB) const;
    /**
     * @brief Calculates the gamma (electron repulsion) matrix.
     */
    void calculate_gamma();
    /**
     * @brief Builds the core Hamiltonian matrix for the molecule.
     */
    void build_core_hamiltonian();
    /**
     * @brief Returns the number of basis functions in the molecule.
     */
    int getN() const;
    /**
     * @brief Returns the number of valence electrons in the molecule.
     */
    int getValElectrons() const;

    // Hartree-Fock matrix functions (hartreeFock)
    /**
     * @brief Returns the diagonal value of the Fock matrix for a basis function.
     */
    double diagonal_value(int mu, bool isAlpha) const;
    /**
     * @brief Returns the off-diagonal value of the Fock matrix for two basis functions.
     */
    double offDiagonal_value(int mu, int nu, bool isAlpha) const;
    /**
     * @brief Builds the Fock matrix for the molecule based on the electron density.
     */
    arma::mat build_fock_matrix(bool isAlpha);
    /**
     * @brief Calculates the number of alpha and beta electrons based on the molecule's configuration.
     */
    void calculate_alpha_beta_electrons();
    /**
     * @brief Calculates the two-electron integral for two primitive Gaussians.
     */
    double I2e_pG(arma::vec &Ra, arma::vec &Rb, double sigmaA, double sigmaB) const;
    /**
     * @brief Evaluates the two-electron integral for atomic orbitals of type s.
     */
    double eval_gamma(const BasisFunction& basisFuncA, const BasisFunction& basisFuncB) const;
    /**
     * @brief Calculates the total electron density for each atom in the molecule.
     */
    void calculate_pTot();
    /**
     * @brief Prints the SCF results for each iteration during the Hartree-Fock calculation.
     */
    void printSCF(bool converged, int it) const;
};

void calculateBondEnergy(std::unordered_map<std::string, double>& energies);

std::string getMoleculeName(const std::string& filePath);

#endif