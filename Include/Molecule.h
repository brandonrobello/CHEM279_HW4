// Chem279: Numerical Algorithms Applied to Computational Quantum Chemistry  
// Creator: Brandon Robello
// Date Created: 10/18/24
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
    double IA_s;  // Value for 1/2 (Is + As)
    double IA_p;  // Value for 1/2 (Ip + Ap)
    double beta;   // Value for -β
};


/**
 * @class Molecule
 * 
 * @brief Represents a Molecule object with its vector of Atoms and associated basis functions.
 * 
 * The Molecule class models a molecular system in quantum chemistry, storing information about its
 * atoms, basis functions, and quantum matrices such as the overlap matrix, Hamiltonian matrix, and molecular
 * orbital coefficients. It provides functionality to compute overlaps between Gaussian basis functions,
 * build the Hamiltonian, and calculate the total energy of the molecule.
 * 
 * The class also allows for the construction of atomic and basis function information,
 * with methods to print the overlap and Hamiltonian matrices.
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
     * @brief Returns the number of basis functions in the molecule.
     * @return int Number of basis functions.
     */
    int getN() const;

    /**
     * @brief Returns the number of valence electrons in the molecule.
     * @return int Number of valence electrons.
     */
    int getValElectrons() const;

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
     * @brief Prints the 
     */
    void printGamma() const;

    /**
     * @brief Prints the orthogonalized overlap matrix (X) matrix of the molecule.
     */
    void printX() const;

    /**
     * @brief Prints the Molecular Orbital Coefficients matrix (C) matrix of the molecule.
     */
    void printC() const;

    /**
     * @brief Prints the Molecular Orbital Overlap matrix (MO_S) of the molecule.
     */ 
    void printMO_S() const;

        /**
     * @brief Print
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
     * 
     * This method computes the full overlap matrix for the molecule, storing it in S_overlap_matrix_.
     */
    void calculate_overlap_matrix();

    /**
     * @brief Builds the Hamiltonian matrix for the molecule.
     * 
     * Constructs the Hamiltonian matrix for the molecule and stores it in H_matrix_.
     */
    void build_hamiltionian();

    /**
     * @brief Constructs the orthogonalized X matrix for the molecule.
     * 
     * This method computes the X matrix, which is used to orthogonalize the overlap matrix (S_overlap_matrix_).
     * The X matrix is essential for transforming the Hamiltonian into an orthogonal basis, which simplifies
     * solving the eigenvalue problem.
     */
    void build_X();

    /**
     * @brief Constructs the molecular orbital coefficient matrix (C_matrix_).
     * 
     * This method computes the molecular orbital coefficients (C_matrix_) by solving the
     * eigenvalue equation for the Hamiltonian matrix in the orthogonal basis. The resulting
     * coefficients describe the molecular orbitals as linear combinations of the basis functions.
     */
    void build_C();

/**
     * @brief Constructs the molecular orbital overlap matrix (MO_S_matrix_).
     * 
     * This method calculates the overlap matrix in the molecular orbital basis (MO_S_matrix_) 
     * by transforming the overlap matrix from the atomic orbital basis using the molecular 
     * orbital coefficients.
     */
    void build_MO_S();

    /**
     * @brief Calculates the total electronic energy of the molecule.
     * 
     * This function computes the total electronic energy of the molecule using the Hamiltonian and
     * the molecular orbital coefficients. The total energy includes the kinetic energy of the electrons
     * and their interactions with the nuclei and each other.
     * 
     * @return double The total electronic energy of the molecule.
     */
    double calc_totalE();

    void build_Fmatrix_wSCF();

    void calculate_gamma();

    void build_core_hamiltonian();

private:
    /**
     * @brief Builds the basis functions for a given atom.
     * 
     * This private method initializes the basis functions for an atom and adds them to the molecule's basis function list.
     * 
     * @param atom The Atom object for which to build basis functions.
     */
    void buildBasisFunctions();
    
    /**
     * @brief Initializes molecule properties, including valence electrons and basis functions.
     * 
     * Called by the constructor to set up the initial state of the molecule.
     */
    void Molecule_init();                       

    std::vector<Atom> atoms_;                   ///< Vector containing the atoms of the molecule.
    std::vector<BasisFunction> basisFunctions_; ///< Store the list basis functions
    arma::mat S_overlap_matrix_;                ///< Overlap Matrix
    arma::mat FA_matrix_;                       ///< Alpha Fock Matrix
    arma::mat FB_matrix_;                       ///< Beta Fock Matrix
    arma::mat pA_matrix_;                       ///< Alpha Denisty Matrix
    arma::mat pB_matrix_;                       ///< Beta Denisty Matrix
    arma::mat H_core_matrix_;
    arma::mat pTot_mu_nu_;
    arma::vec pTot_;                            ///< Total density vector for each atom
    arma::mat VA_matrix_;                       ///< Alpha Eigen vectors Matrix
    arma::mat VB_matrix_;                       ///< Beta Eigen vectors Matrix
    arma::vec epsilon_A_;                        ///< Alpha Eigen values vector
    arma::vec epsilon_B_;                        ///< Beta Eigen values vector
    arma::mat cA_matrix_;                       ///< Alpha MO coefficients Matrix
    arma::mat cB_matrix_;                       ///< Beta MO coefficients Matrix
    arma::mat gamma_matrix_;
    int N_ = 0;                                 ///< Number of basis functions in the molecule.
    int val_electrons_ = 0;                     ///< Number of valence electrons in the molecule
    int numCarbons_ = 0;                        ///< Number of carbon atoms
    int numHydrogens_ = 0;                      ///< Number of hydrogen atoms
    int alpha_e_, beta_e_ = 0;                  ///< num of p and q electrons in the MO
    int alpha_electron_E_, beta_electron_E_ = 0;                  ///< 
    int charge_ = 0;                            ///< Charge of molecule
    double nuclear_repulsion_E_ = 0;
    std::unordered_map<int, CNDO2Parameters> semiEmpiricalParams_; ///< Map with semi-empirical parameters

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
    

    // Returns the diagonal value of the fock matrix for a basis function and current electron density.
    double diagonal_value(int mu, bool isAlpha) const;

    // Returns the off-diagonal value of the fock matrix for two basis functions.
    double offDiagonal_value(int mu, int nu, bool isAlpha) const;

    // Builds the Fock matrix for the molecule.
    arma::mat build_fock_matrix(bool isAlpha);

    /**
     * @brief Get the p and q numbers for the valence electrons and charge of MO
     */
    void calculate_alpha_beta_electrons();

    // 2-electron integral of two primitive Gaussians
    double I2e_pG(arma::vec &Ra, arma::vec &Rb, double sigmaA, double sigmaB) const;
    // Function to evaluate 2-electron integral for atomic orbitals (s-type)
    double eval_gamma(const BasisFunction& basisFuncA, const BasisFunction& basisFuncB) const;

    void calculate_pTot();

    void print_armaMat(const std::string name, const arma::mat& mat) const;
};

#endif