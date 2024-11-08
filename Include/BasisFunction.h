// Chem279: Numerical Algorithms Applied to Computational Quantum Chemistry  
// Creator: Brandon Robello
// Date Created: 10/18/24
//
// BasisFunction.h contains the necessary include files and 
//        C++ Class declarations for BasisFunction class

#include "Atom.h"
#include <armadillo>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
#include <set>
#include <map>

#ifndef BASISFUNC_H
#define BASISFUNC_H

/**
 * @class BasisFunction
 * 
 * @brief Represents a contracted Gaussian-type orbital (GTO) as a basis function for atoms in a molecule.
 * 
 * This class is responsible for modeling the basis function used in quantum chemistry computations,
 * particularly in molecular orbital theory. It stores the properties of a Gaussian function for a given
 * atom, such as the atomic number, center coordinates, angular momenta, exponents (alphas),
 * contraction coefficients, and normalization constants. These properties are essential in the
 * evaluation of atomic integrals for quantum chemical calculations.
 * 
 * The basis function is defined by a set of Gaussian primitives that are combined into a single
 * contracted function. The normalization constants for each primitive are precomputed based on
 * the angular momenta and exponents provided.
 * 
 * @note Each basis function is associated with a specific atom in the molecule, and the basis functions
 * for different atoms are combined to build the molecular wavefunctions.
 * 
 * @details 
 * - **element**: Represents the atomic number of the atom associated with the basis function.
 * - **center**: The 3D Cartesian coordinates (X, Y, Z) of the atom at which the basis function is centered.
 * - **ang_momenta**: The angular momentum quantum numbers (l, m, n) that define the shape of the Gaussian.
 * - **alphas**: A vector of exponents that determine the width of each Gaussian primitive in the contraction.
 * - **coeffs**: The contraction coefficients that weight each Gaussian primitive in the overall basis function.
 * - **norms**: Precomputed normalization constants for each Gaussian primitive, ensuring the function is normalized.
 * 
 * @see Molecule, Atom
 */
class BasisFunction {
public:
    int element;                                ///< Atomic Number
    int val_electrons;
    int orbitalType;
    int atomic_pos;
    arma::vec center;                           ///< Center (X, Y, Z) of the atom
    arma::vec ang_momenta;                      ///< Angular momenta
    std::vector<double> alphas;                 ///< Exponents
    std::vector<double> coeffs;                 ///< Contraction coefficients
    std::vector<double> norms;                  ///< Normalization constants

    /**
     * @brief Constructor for BasisFunction class.
     * 
     * This constructor initializes a basis function with its corresponding
     * atomic number, center coordinates, angular momenta, exponents (alphas),
     * and contraction coefficients. It also calculates the normalization constants
     * for each exponent based on the angular momenta.
     * 
     * @param E Atomic number of the element
     * @param R Center coordinates (X, Y, Z) of the atom
     * @param L Angular momenta (l, m, n)
     * @param a Exponents (alphas) for the Gaussian functions
     * @param d Contraction coefficients (coeffs)
     * @param pos
     */
    BasisFunction(const int pos, const int E, const arma::vec& R, arma::vec L, 
                    const std::vector<double>& a, const std::vector<double>& d)
        : atomic_pos(pos), element(E), center(R), ang_momenta(L), alphas(a), coeffs(d) {
        // Compute normalization constants
        for (size_t i = 0; i < alphas.size(); ++i) {
            double norm = computeNormalization(alphas[i], ang_momenta);
            norms.push_back(norm);
        }

        // Set Orbital Type
        orbitalType = arma::accu(ang_momenta);

        // Note Val_electrons
        if (E == 6) { // Carbon contributes 4 valence electrons
            val_electrons = 4;
        } else if (E == 1) { // Hydrogen contributes 1 electron
            val_electrons = 1;
        } else if (E == 7) { // Nitrogen contributes 4 orbital basis functions  and 5 electrons
            val_electrons = 5;
        } else if (E == 8) { // Oxygen contributes 4 orbital basis functions  and 6 electrons
            val_electrons = 6;
        } else if (E == 9) { // Fluorine contributes 4 orbital basis functions  and 7 electrons
            val_electrons = 7;
        } else {
            throw std::invalid_argument("Unsupported element with atomic number: " + std::to_string(element));
        }
    }

    /**
     * @brief Get the angular momentum vector (l, m, n) for this basis function.
     * 
     * This function returns the angular momenta of the basis function, which 
     * defines the type of orbital represented (e.g., s, p, d orbitals in quantum chemistry).
     * 
     * @return arma::vec Angular momentum quantum numbers (l, m, n).
     */
    arma::vec getL() const {return ang_momenta;}

    /**
     * @brief Overload << operator for printing BasisFunction details
     * @param os The output stream
     * @param bf The BasisFunction to be printed
     * @return The output stream with the BasisFunction's details
     */
    friend std::ostream& operator<<(std::ostream& os, const BasisFunction& bf) {
    os << "Basis Function Center: (" << bf.center(0) << ", " << bf.center(1) << ", " << bf.center(2) << ")\n";
    os << "Angular Momentum (l, m, n): (" << bf.ang_momenta(0) << ", " << bf.ang_momenta(1) << ", " << bf.ang_momenta(2) << ")\n";
    os << "Exponents (alphas): ";
    for (const auto& alpha : bf.alphas) {
        os << alpha << " ";
    }
    os << "\n";
    os << "Contraction Coefficients (coeffs): ";
    for (const auto& coeff : bf.coeffs) {
        os << coeff << " ";
    }
    os << "\n";
    os << "Normalization Constants (norms): ";
    for (const auto& norm : bf.norms) {
        os << norm << " ";
    }
    os << "\n";
    return os;
    }
private:
    /**
     * @brief Compute the normalization constant for a Gaussian primitive.
     * 
     * This function calculates the normalization constant for a Gaussian primitive,
     * based on its exponent and angular momentum quantum numbers. The normalization ensures
     * that the Gaussian-type orbital is properly scaled in quantum chemical calculations.
     * 
     * @param alpha Exponent of the Gaussian function.
     * @param ang_momenta Angular momentum vector (l, m, n) for the basis function.
     * @return double The normalization constant for the Gaussian primitive.
     */
    double computeNormalization(double alpha, arma::vec ang_momenta) {
        return pow(2 * alpha / M_PI, 3.0 / 4.0) *
               pow(4 * alpha, arma::accu(ang_momenta) / 2.0) / 
               (factorial(2*ang_momenta(0) - 1) * factorial(2*ang_momenta(1) - 1) * factorial(2*ang_momenta(2) - 1));
    }

    /**
     * @brief Compute the double factorial (n!!) of an integer.
     * 
     * This helper function calculates the double factorial of a number, 
     * which is useful in normalization constants for angular momentum. 
     * For example, the double factorial of 5 is 5 * 3 * 1 = 15.
     * 
     * @param n Integer input to compute the double factorial.
     * @return int The result of the double factorial computation.
     */
    int factorial(int n) {
        return (n <= 1) ? 1 : n * factorial(n - 2);
    }
};

#endif