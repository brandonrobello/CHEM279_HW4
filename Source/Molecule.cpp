// Chem 279: Numerical Algorithms Applied to Computational Quantum Chemistry            
// Creator: Brandon Robello
// Date Created: 10/18/24
//    
// This file contains C++ implementation of a Molecule. 
// The API (include) file Molecule.h, which includes 
// prototypes of the definitions in this file.


#include "Atom.h"
#include "Molecule.h"
#include "BasisFunction.h"

#include <cmath>
#include <stdexcept>
#include <cassert>
#include <iostream>
#include <armadillo>

#include <unordered_map>
#include <string>

#include <iostream>
#include <vector>
#include <algorithm> // For std::count

// Function to generate electron configuration
std::vector<int> generateElectronConfiguration(int totalElectrons) {
    // Electron capacity of common orbitals (1s, 2s, 2p, 3s, 3p, 4s, 3d, 4p, 5s, etc.)
    std::vector<int> orbitalCapacities = {2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6};
    std::vector<int> electronConfiguration;

    int remainingElectrons = totalElectrons;

    for (int capacity : orbitalCapacities) {
        if (remainingElectrons > 0) {
            int electronsInOrbital = std::min(remainingElectrons, capacity);
            electronConfiguration.push_back(electronsInOrbital);
            remainingElectrons -= electronsInOrbital;
        } else {
            break;
        }
    }

    return electronConfiguration;
}

// Function to count unpaired electrons in a given electron configuration
int countUnpairedElectrons(const std::vector<int>& electronConfiguration) {
    int unpairedElectrons = 0;

    // Loop through each orbital in the configuration
    for (int electrons : electronConfiguration) {
        // Check if the orbital has unpaired electrons
        if (electrons % 2 != 0) {
            unpairedElectrons++;
        }
    }

    return unpairedElectrons;
}

// Function to determine the multiplicity based on the number of unpaired electrons
int calculateMultiplicity(int unpairedElectrons) {
    int S = unpairedElectrons / 2; // Total spin quantum number S
    int multiplicity = 2 * S + 1;  // Multiplicity formula: M = 2S + 1

    return multiplicity;
}

// Function to determine the number of alpha and beta electrons based on atoms of the molecule
void Molecule::calculate_alpha_beta_electrons() {
    int totalElectrons = 0;

    for (auto& atom : atoms_){
        totalElectrons += atom.getValElectrons();
    }

    std::vector<int> electronConfiguration = generateElectronConfiguration(totalElectrons);
    int unpairedElectrons = countUnpairedElectrons(electronConfiguration);
    int multiplicity = calculateMultiplicity(unpairedElectrons);


    // Calculate total spin quantum number S from multiplicity (M = 2S + 1 -> S = (M - 1) / 2)
    int spinQuantumNumber = (multiplicity - 1) / 2;

    // Calculate the number of alpha and beta electrons
    alpha_e_ = static_cast<int>(std::ceil((totalElectrons + spinQuantumNumber) / 2.0));
    beta_e_ = static_cast<int>(std::floor((totalElectrons - spinQuantumNumber) / 2.0));

    if (alpha_e_ < 0 || beta_e_ < 0) {
        throw std::runtime_error("Calculated alpha or beta electrons are negative, check your input and calculations.");
    }

    std::cout << "Alpha electrons (p): " << alpha_e_ << std::endl;
    std::cout << "Beta electrons (q): " << beta_e_ << std::endl;
}



// Constructor and intializing Functions
Molecule::Molecule(const std::vector<Atom>& Atoms)
    : atoms_(Atoms){
    // Map with semi-empirical parameters for H, C, N, O, F
    semiEmpiricalParams_ = {
        {1, {7.176, 0.0, -9}},    // H: Ionization energy and beta (only s-orbital)
        {6, {14.051, 5.572, -21}}, // C: Ionization and p-orbital values
        {7, {19.316, 7.275, -25}}, // N: Ionization and p-orbital values
        {8, {25.390, 9.111, -31}}, // O: Ionization and p-orbital values
        {9, {32.272, 11.080, -39}} // F: Ionization and p-orbital values
    };

    // Initalize process for MO evaluation
    Molecule_init();
}

// Builds the basis functions for a given atom.
void Molecule::buildBasisFunctions() {
    
    int pos = 0;
    for (auto atom : atoms_){
        int E = atom.getAtomicNumber();
        std::vector<std::vector<double>>  coeffs;
        arma::vec center = { atom.getX(), atom.getY(), atom.getZ() };
        arma::vec L_s = { 0, 0, 0 };    // s: l=m=n=0
        arma::vec L_px = { 1, 0, 0 };   // px: l=1, m=n=0
        arma::vec L_py = { 0, 1, 0 };   // py: l=0, m=1, n=0
        arma::vec L_pz = { 0, 0, 1 };   // pz: l=m=0, n=1

        if (E == 1) {
            coeffs = {
                { 0.15432897, 0.53532814, 0.44463454 } // 1s orbital
            };
        }
        else {
            coeffs = {
                { -0.09996723, 0.39951283, 0.70011547 }, // 2s orbital
                { 0.15591627, 0.60768372, 0.39195739 }   // 2p orbitals
            };
        }

        if (E == 1) {
            // Basis functions for H (1s orbital)
            std::vector<double> alphas = { 3.42525091, 0.62391373, 0.16885540 };
            basisFunctions_.emplace_back(pos, E, center, L_s, alphas, coeffs[0]);
        } else if (E == 6) {
            // Basis functions for C (2s and 2p orbitals)
            std::vector<double> alphas = { 2.94124940, 0.68348310, 0.22228990 };

            basisFunctions_.emplace_back(pos, E, center, L_s, alphas, coeffs[0]);
            basisFunctions_.emplace_back(pos, E, center, L_px, alphas, coeffs[1]);
            basisFunctions_.emplace_back(pos, E, center, L_py, alphas, coeffs[1]);
            basisFunctions_.emplace_back(pos, E, center, L_pz, alphas, coeffs[1]);
        } else if (E == 7) {
            // Basis functions for N (2s and 2p orbitals)
            std::vector<double> alphas = { 3.78045590, 0.87849660, 0.28571440 };

            basisFunctions_.emplace_back(pos, E, center, L_s, alphas, coeffs[0]);
            basisFunctions_.emplace_back(pos, E, center, L_px, alphas, coeffs[1]);
            basisFunctions_.emplace_back(pos, E, center, L_py, alphas, coeffs[1]);
            basisFunctions_.emplace_back(pos, E, center, L_pz, alphas, coeffs[1]);
        } else if (E == 8) {
            // Basis functions for O (2s and 2p orbitals)
            std::vector<double> alphas = { 5.03315130, 1.16959610, 0.38038900 };

            basisFunctions_.emplace_back(pos, E, center, L_s, alphas, coeffs[0]);
            basisFunctions_.emplace_back(pos, E, center, L_px, alphas, coeffs[1]);
            basisFunctions_.emplace_back(pos, E, center, L_py, alphas, coeffs[1]);
            basisFunctions_.emplace_back(pos, E, center, L_pz, alphas, coeffs[1]);
        } else if (E == 9) {
            // Basis functions for F (2s and 2p orbitals)
            std::vector<double> alphas = { 6.46480320, 1.50228120, 0.48858850 };

            basisFunctions_.emplace_back(pos, E, center, L_s, alphas, coeffs[0]);
            basisFunctions_.emplace_back(pos, E, center, L_px, alphas, coeffs[1]);
            basisFunctions_.emplace_back(pos, E, center, L_py, alphas, coeffs[1]);
            basisFunctions_.emplace_back(pos, E, center, L_pz, alphas, coeffs[1]);
        } else {
            // Error handling for unsupported atoms
            throw std::invalid_argument("Error: Basis functions for atomic number " + std::to_string(E) + " are not defined.");
        }

        pos++;
    }
}


// Initializes molecule properties, including valence electrons and basis functions.
void Molecule::Molecule_init(){
    // Get the basis functions and valence electrons for each atom
    for (auto& atom : atoms_) {
        N_ += atom.getN();
        val_electrons_ += atom.getValElectrons();
    }

    // Build Basis Functions
    buildBasisFunctions();

    // Get the p and q numbers for the valence electrons and charge of MO
    calculate_alpha_beta_electrons();

    // // Check the condition: val_electrons should equal 2*C + H/2
    // int expected_pairs = 2 * numCarbons_ + numHydrogens_ / 2;
    // if (val_electrons_/2 != expected_pairs) {
    //     std::cout << "val_electrons_: " << val_electrons_ << "\nexpected_pairs: " << expected_pairs << std::endl;
    //     throw std::runtime_error("Error: Number of valence electrons does not satisfy the condition 2*C + H/2");
    // }

    // Set S_overlap_matrix_ to be a N x N matrix and initialize it to zeros
    S_overlap_matrix_.set_size(N_, N_);
    S_overlap_matrix_.zeros();

    // Evaluate the S_overlap_matrix
    calculate_overlap_matrix();

    // Initialize alpha, beta, and total density matrices to zero (initial guess)
    pA_matrix_.set_size(N_, N_);
    pB_matrix_.set_size(N_, N_);
    pTot_.set_size(atoms_.size());
    cA_matrix_.set_size(N_, N_);
    cB_matrix_.set_size(N_, N_);
    pA_matrix_.zeros();
    pB_matrix_.zeros();
    pTot_.zeros();
    cA_matrix_.zeros();
    cB_matrix_.zeros();

    // Initialize alpha and beta fock matrices to zero
    FA_matrix_.set_size(N_, N_);
    FB_matrix_.set_size(N_, N_);
    FA_matrix_.zeros();
    FB_matrix_.zeros();

    // Build Fock matrix with Unrestricted Hartree-Fock CNDO/2 method
    // build_Fmatrix_wSCF();
}

// Adds a new atom to the molecule.
void Molecule::addAtom(const Atom& atom) {
    atoms_.push_back(atom);

    // Reset the private members
    N_ = 0;
    val_electrons_ = 0;
    numCarbons_ = 0;
    numHydrogens_ = 0;

    // Clear the basisFunctions_ vector to remove the previous basis functions
    basisFunctions_.clear();

    // Reinitialize the molecule
    Molecule_init();
}

// Get Functions
int Molecule::getN() const {return N_;}
int Molecule::getValElectrons() const {return val_electrons_;}
std::vector<Atom> Molecule::getAtoms() const {return atoms_;}

// Print Functions
std::ostream& operator<<(std::ostream& os, const Molecule& molecule) {
    os << "Basis Functions (N): " << molecule.getN() 
    <<  "\nValence Elections (2n): " << molecule.getValElectrons() 
    <<  "\nAtoms:\n";
    for (auto& atom : molecule.getAtoms())
    os << atom << "\n";
    return os;
}

void Molecule::printBasisFunctions() const {
    int iter = 0;
    for (auto& func : basisFunctions_){
        std::cout << "BasisFunctions " << iter << ":\n" << func << std::endl;
        iter++;
    }
}

void Molecule::printS() const {
    std::cout << "Overlap Matrix:" << std::endl;
    std::cout << S_overlap_matrix_ << std::endl;
}

void Molecule::printH() const {
    std::cout << "H_core Matrix:" << std::endl;
    std::cout << H_core_matrix_ << std::endl;
}

void Molecule::printGamma() const {
    std::cout << "Gamma Matrix:" << std::endl;
    std::cout << gamma_matrix_ << std::endl;
}

void Molecule::print_armaMat(const std::string name, const arma::mat& mat) const {
    std::cout << name << " Matrix:" << std::endl;
    std::cout << mat << std::endl;
}

// Function to calculate the exponential prefactor in Eq. 2.9
double Molecule::exponential_prefactor(double alpha, double beta, double X_A, double X_B) const {
    double exponent = (-alpha * beta * std::pow((X_A - X_B), 2)) / (alpha + beta);

    double root_term = std::sqrt(M_PI / (alpha + beta));

    return std::exp(exponent) * root_term;
}

// Function to calculate factorial (n!)
int Molecule::factorial(int n) const {
    if (n < 0) {
        throw std::invalid_argument("Factorial is not defined for negative numbers.");
    }
    int fact = 1;
    for (int i = 2; i <= n; ++i) {
        fact *= i;
    }
    return fact;
}

// Function to calculate double factorial (n!!)
int Molecule::double_factorial(int n) const {
    if (n < -1) {
        throw std::invalid_argument("Double factorial is not defined for negative numbers.");
    }
    if (n == 0 || n == -1) return 1;  // Handle 0!! = 1 and -1!! = 1 by definition
    int fact = 1;
    for (int i = n; i > 0; i -= 2) {
        fact *= i;
    }
    return fact;
}

// Function to calculate binomial coefficient
double Molecule::calculate_binomial(int l, int i) const {
    if(l < i || i < 0){
        throw std::invalid_argument("Comination: the number of elements should be bigger than selection numbers AND two numbers should be positive\n");
        return EXIT_FAILURE;
    }

    return factorial(l) / (factorial(i) * factorial(l - i));
}

// Function to calculate the center of the product Gaussian
double Molecule::calculate_RP(double alpha, double beta, double X_A, double X_B) const {
    return ((alpha * X_A) + (beta * X_B)) / (alpha + beta);
}

// Helper function to calculate the power terms with edge case handling
double Molecule::calculate_power_term(double Rp, double center, int l, int i) const {
    return std::pow(Rp - center, l - i);
}

// Function to calculate the term for in eq. 2.9 for the overlap matrix in one dimension calculation
double Molecule::calculate_term(double alpha, double beta, double X_A, double X_B, int lA, int lB, int i, int j) const {
    int double_fact = double_factorial(i + j - 1);
    double Rp = calculate_RP(alpha, beta, X_A, X_B);

    // Calculate the power terms using the helper function
    double pow_XA = calculate_power_term(Rp, X_A, lA, i);
    double pow_XB = calculate_power_term(Rp, X_B, lB, j);

    // Compute the numerator
    double num = double_fact * pow_XA * pow_XB;

    // Compute the denominator
    double denom = std::pow(2 * (alpha + beta), double(i + j) / 2.0);

    return num / denom;
}

// Calculates the one-dimensional overlap integral between two Gaussian primitives.
double Molecule::Overlap_onedim(double alpha, double beta, double X_A, double X_B,  int lA, int lB) const
{   // Exponential prefactor outside the summation
    double prefactor = exponential_prefactor(alpha, beta, X_A, X_B);

    double result = 0.0;
    for(int i = 0; i<= lA; i++)
        for(int j = 0; j <= lB; j++){
        // Skip odd (i + j) terms
        if(( i + j) % 2 == 1)
            continue;
        
        // Calculate the binomial coefficients (lA choose i) and (lB choose j)
        int binomial_A = calculate_binomial(lA, i);
        int binomial_B = calculate_binomial(lB, j);

        // Calculate the term in the summation
        double term = calculate_term(alpha, beta, X_A, X_B, lA, lB, i, j);

        double temp = binomial_A * binomial_B * term;

        result += temp;
        }

    result *= prefactor;
    return result;
}

// Calculates the primitive overlap between two Gaussian primitives.
double Molecule::calculate_primitive_overlap(double alpha_k, double alpha_l, const arma::vec& R_A, const arma::vec& R_B, const arma::vec lA, const arma::vec lB) const {
    // Calculate overlap in x, y, and z dimensions
    double S_x = Overlap_onedim( alpha_k, alpha_l, R_A(0), R_B(0), lA(0), lB(0));  // x-direction
    double S_y = Overlap_onedim( alpha_k, alpha_l, R_A(1), R_B(1), lA(1), lB(1));  // y-direction
    double S_z = Overlap_onedim( alpha_k, alpha_l, R_A(2), R_B(2), lA(2), lB(2));  // z-direction

    // Return total overlap as the product of the one-dimensional overlaps
    return S_x * S_y * S_z;
}

// Calculates the contracted overlap between two basis functions.
double Molecule::calculate_contracted_overlap(int mu, int nu) const {
    double S_mu_nu = 0.0;
    
    // Loop over primitives in basis function mu
    for (int k = 0; k < 3; ++k) {
        double d_mu_k = basisFunctions_[mu].coeffs[k];
        double N_mu_k = basisFunctions_[mu].norms[k];
        double alpha_mu_k = basisFunctions_[mu].alphas[k];

        // Loop over primitives in basis function nu
        for (int l = 0; l < 3; ++l) {

            double d_nu_l = basisFunctions_[nu].coeffs[l];
            double N_nu_l = basisFunctions_[nu].norms[l];
            double alpha_nu_l = basisFunctions_[nu].alphas[l];

            // Get the centers of the primitives
            arma::vec R_A = basisFunctions_[mu].center;  // Coordinates of center of primitive k
            arma::vec R_B = basisFunctions_[nu].center;  // Coordinates of center of primitive l

            // Get the angular momentum of the primitives
            arma::vec lA = basisFunctions_[mu].getL();
            arma::vec lB = basisFunctions_[nu].getL();

            // Calculate the primitive overlap integral S_kl
            double S_kl = calculate_primitive_overlap(alpha_mu_k, alpha_nu_l, R_A, R_B, lA, lB);


            // Add the contribution to the total overlap
            S_mu_nu += d_mu_k * d_nu_l * N_mu_k * N_nu_l * S_kl;
        }
    }

    return S_mu_nu;  // Return the total overlap integral for the contracted basis functions
}

// Function to calculate overlay matrix
void Molecule::calculate_overlap_matrix() {
    // Loop over all pairs of basis functions (mu, nu)
    for (int mu = 0; mu < N_; ++mu) {
        for (int nu = 0; nu < N_; ++nu) {
            // Calculate the contracted overlap integral between basis function mu and nu
            double S_mu_nu = calculate_contracted_overlap(mu, nu);

            // Store the result in the overlap matrix
            S_overlap_matrix_(mu, nu) = S_mu_nu;
        }
    }
}

// 2-electron integral of two primitive Gaussians
double Molecule::I2e_pG(arma::vec &Ra, arma::vec &Rb, double sigmaA, double sigmaB) const {
    double conversion_factor = 27.211324570273; // a.u. to eV

    double V2 = 1.0 / (sigmaA + sigmaB); // Equation 3.9
    double Rd = arma::norm(Ra - Rb, 2);
    double T = V2 * std::pow(Rd, 2); // Equation 3.7

    double UA = std::pow(M_PI * sigmaA, 1.5); // Equation 3.11
    double UB = std::pow(M_PI * sigmaB, 1.5); // Equation 3.11

    if (Rd == 0.0) {
        // If RA == RB, use equation 3.15
        double result = UA * UB * 2.0 * std::sqrt(V2 / M_PI);
        return result * conversion_factor; // Converting a.u. to eV consistent with the semi-empirical params
    }

    // Equation 3.14
    double result = UA * UB * std::sqrt(1.0 / (Rd * Rd)) * std::erf(std::sqrt(T));

    return result * conversion_factor; // Converting a.u. to eV consistent with the semi-empirical params
}

// Function to evaluate 2-electron integral for atomic orbitals (s-type)
double Molecule::eval_gamma(const BasisFunction& basisFuncA, const BasisFunction& basisFuncB) const {
    int len = basisFuncA.alphas.size();
    assert(basisFuncB.alphas.size() == len);

    arma::vec Ra = basisFuncA.center;
    arma::vec Rb = basisFuncB.center;
    const std::vector<double>& alpha_a = basisFuncA.alphas;
    const std::vector<double>& da = basisFuncA.coeffs;
    const std::vector<double>& norm_a = basisFuncA.norms;
    const std::vector<double>& alpha_b = basisFuncB.alphas;
    const std::vector<double>& db = basisFuncB.coeffs;
    const std::vector<double>& norm_b = basisFuncB.norms;


    double gamma = 0.0;
    // Loop over k, k', l, l' as in Equation 3.3
    for (size_t k1 = 0; k1 < len; ++k1){
        for (size_t k2 = 0; k2 < len; ++k2) {
            double sigmaA = 1.0 / (alpha_a[k1] + alpha_a[k2]); // Equation 3.10
            for (size_t l1 = 0; l1 < len; ++l1) {
                for (size_t l2 = 0; l2 < len; ++l2) {
                    double sigmaB = 1.0 / (alpha_b[l1] + alpha_b[l2]); // Equation 3.10
                    double I2e = I2e_pG(Ra, Rb, sigmaA, sigmaB); // Integral calculation
                    gamma += (da[k1] * norm_a[k1]) * (da[k2] * norm_a[k2]) * 
                            (db[l1] * norm_b[l1]) * (db[l2] * norm_b[l2]) * I2e; // Equation 3.3
                }
            }
        }
    }
    return gamma;
}

// Returns the diagonal value of the fock matrix for a basis function and current electron density.
double Molecule::diagonal_value(int mu, bool isAlpha) const {
    // Set the density matrix
    const arma::mat& densityMatrix = isAlpha ? pA_matrix_ : pB_matrix_;

    // Get element atomic position
    int element = basisFunctions_[mu].element;
    int A = basisFunctions_[mu].atomic_pos;

    if (semiEmpiricalParams_.find(element) == semiEmpiricalParams_.end()) {
        throw std::invalid_argument("Element not supported in semi-empirical parameter map.");
    }

    // Get semiempirical values
    CNDO2Parameters params = semiEmpiricalParams_.at(element);

    // Determine the correct ionization/affinity energies based on orbital type
    double IA_energy;
    if (basisFunctions_[mu].orbitalType == 0 ) { // S-orbital
        IA_energy = semiEmpiricalParams_.at(element).IA_s;
    } else if (basisFunctions_[mu].orbitalType == 1) { // P-orbital
        IA_energy = semiEmpiricalParams_.at(element).IA_p;
    } else {
        throw std::runtime_error("Unknown orbital type in basis function.");
    }

    // // Debug
    // std::cout << "A = " << A << std::endl;

    // Compute the total electron density minus the valence electrons
    double ptot_Za = pTot_(A) - atoms_[A].getValElectrons();

    // // Debug
    // std::cout << "atoms_[A].getValElectrons() = " << atoms_[A].getValElectrons() << std::endl;

    // Compute either p^alpha or p^beta - 1/2
    double pA_uu_half = densityMatrix(mu, mu) - 0.5;

    // Calculate the self-repulsion term (AA)
    double gamma_AA = gamma_matrix_(A, A);
    double eRepulsion_AA_exchange = (ptot_Za - pA_uu_half) * gamma_AA;

    // Loop to calculate gamma_AB terms for interactions with other atoms
    double sum_gamma_AB = 0.0;
    for (int B = 0; B < atoms_.size(); ++B) {
        if (B != A) {
            sum_gamma_AB += (pTot_(B) - atoms_[B].getValElectrons()) * gamma_matrix_(A, B);
        }
    }

    // Equation 1.4 to calculate the final diagonal Fock matrix element
    return -IA_energy + eRepulsion_AA_exchange + sum_gamma_AB;
}


// Returns the off-diagonal value of the fock matrix for two basis functions.
double Molecule::offDiagonal_value(int mu, int nu, bool isAlpha) const{
// NEED TO EDIT USING Eq 1.5
    // Set the density matrix
    const arma::mat& densityMatrix = isAlpha ? pA_matrix_ : pB_matrix_;

    // 
    int elementA = basisFunctions_[mu].element;
    int A = basisFunctions_[mu].atomic_pos;
    int elementB = basisFunctions_[nu].element;
    int B = basisFunctions_[nu].atomic_pos;

    // Get semiempirical parameters
    CNDO2Parameters params_A = semiEmpiricalParams_.at(elementA);
    CNDO2Parameters params_B = semiEmpiricalParams_.at(elementB);

    // Get overlap matrix element
    double s_uv = S_overlap_matrix_(mu, nu);

    // Calculate average B
    double beta_sum = params_A.beta + params_B.beta;

    // Calculate p(mu, nu) and gamma_AB
    double p_mu_nu = densityMatrix(mu, nu);
    double gamma_AB = gamma_matrix_(A, B);

    return (0.5 * beta_sum * s_uv) - (p_mu_nu * gamma_AB);
}

// Builds the Fock matrix for the molecule.
arma::mat Molecule::build_fock_matrix(bool isAlpha) {
    arma::mat Fock_matrix;
    Fock_matrix.set_size(N_, N_);

    // // Debug
    // std::cout << "Inside build_fock_matrix: " << std::endl;
    // std::cout << "N_: " << N_ << std::endl;

    for (int mu = 0; mu < N_; ++mu) {
        for (int nu = 0; nu < N_; ++nu) {
            // // Debug
            // std::cout << "mu: " << mu << ", nu: " << nu << std::endl;
        
            // Get Overlap Matrix value
            double S_uv = S_overlap_matrix_(mu, nu);

            // // Debug
            // std::cout << "S_uv = " << S_uv << std::endl;

            // Use parameters for diagonal and off-diagonal calculations as needed
            double h = (mu == nu) ? diagonal_value(mu, isAlpha) : offDiagonal_value(mu, nu, isAlpha);

            // // Debug
            // std::cout << "h = " << h << std::endl;

            // Store in the Fock matrix
            Fock_matrix(mu, nu) = h;
        }
    }
    return Fock_matrix;
}

void Molecule::build_Fmatrix_wSCF() {
    bool converged = false;
    int max_iterations = 100;
    double convergence_threshold = 1e-6;

    // Temporary Matrices
    arma::mat Falpha(N_, N_, arma::fill::zeros);
    arma::mat Fbeta(N_, N_, arma::fill::zeros);

    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        std::cout << "Iteration: " << iteration << std::endl;
        std::cout << "p: " << alpha_e_ << ", q: " << beta_e_ << std::endl;


        bool alpha = true;
        bool beta = false;
        // 1. Build the Fock matrices F_alpha and F_beta based on current density
        Falpha = build_fock_matrix(alpha);
        Fbeta = build_fock_matrix(beta);

        std::cout << "F_alpha Matrix:\n" << Falpha << std::endl;
        std::cout << "F_beta Matrix:\n" << Fbeta << std::endl;

        // 2. Solve F C = S C epsilon for alpha and beta electrons
        arma::eig_sym(epsilon_A_, cA_matrix_, Falpha);
        arma::eig_sym(epsilon_B_, cB_matrix_, Fbeta);

    
        std::cout << "C_alpha Matrix:\n" << cA_matrix_ << std::endl;
        std::cout << "C_beta Matrix:\n" << cB_matrix_ << std::endl;

        // 3. Update the alpha and beta density matrices
        arma::mat new_pA_matrix = cA_matrix_.cols(0, alpha_e_ - 1) * cA_matrix_.cols(0, alpha_e_ - 1).t();
        arma::mat new_pB_matrix = cB_matrix_.cols(0, beta_e_ - 1) * cB_matrix_.cols(0, beta_e_ - 1).t();

        std::cout << "New P_alpha Matrix:\n" << new_pA_matrix << std::endl;
        std::cout << "New P_beta Matrix:\n" << new_pB_matrix << std::endl;

        // 4. Check for convergence
        double diff_alpha = arma::accu(arma::abs(new_pA_matrix - pA_matrix_));
        double diff_beta = arma::accu(arma::abs(new_pB_matrix - pB_matrix_));

        std::cout << "Difference (Alpha): " << diff_alpha << std::endl;
        std::cout << "Difference (Beta): " << diff_beta << std::endl;

        if (diff_alpha < convergence_threshold && diff_beta < convergence_threshold) {
            converged = true;
            break;
        }

        // Update the density matrices for the next iteration
        pA_matrix_ = new_pA_matrix;
        pB_matrix_ = new_pB_matrix;
        pTot_mu_nu_ = pA_matrix_ + pB_matrix_;
        calculate_pTot();

        std::cout << "P_total:\n" << pTot_ << std::endl;

    }   

    if (converged) {
        std::cout << "Eigenvalues (Alpha):\n" << epsilon_A_ << std::endl;
        std::cout << "Eigenvalues (Beta):\n" << epsilon_B_ << std::endl; 

        std::cout << "SCF converged successfully." << std::endl;
        FA_matrix_ = Falpha;
        FB_matrix_ = Fbeta;

    } else {
        std::cerr << "SCF did not converge within the maximum number of iterations." << std::endl;
    }
}

// Function to calculate the gamma matrix for the molecule
void Molecule::calculate_gamma() {
    int int_gamma_size = atoms_.size();
    gamma_matrix_.set_size(int_gamma_size, int_gamma_size);

    // Loop over each pair of basis functions to calculate gamma
    for (int A = 0; A < int_gamma_size; ++A) {
        for (int B = 0; B < int_gamma_size; ++B) {
            gamma_matrix_(A, B) = eval_gamma(basisFunctions_[A], basisFunctions_[B]);
        }
    }
}


// Function to build the core Hamiltonian matrix (H_core)
void Molecule::build_core_hamiltonian() {
    H_core_matrix_.set_size(N_, N_);

    for (int mu = 0; mu < N_; ++mu) {
        for (int nu = 0; nu < N_; ++nu) {
            if (mu == nu) {
                // Diagonal elements (Equation 2.6)
                int A = basisFunctions_[mu].atomic_pos;
                // Get Atomic element
                int element_A = atoms_[A].getAtomicNumber();
            // std::cout << "Element_mu";
                int Z_A = atoms_[A].getValElectrons();
                double gamma_AA = gamma_matrix_(A, A);

                // Determine the correct ionization/affinity energies based on orbital type
                double IA_energy;
                if (basisFunctions_[mu].orbitalType == 0 ) { // S-orbital
                    IA_energy = semiEmpiricalParams_.at(element_A).IA_s;
                } else if (basisFunctions_[mu].orbitalType == 1) { // P-orbital
                    IA_energy = semiEmpiricalParams_.at(element_A).IA_p;
                } else {
                    throw std::runtime_error("Unknown orbital type in basis function.");
                }
            // std::cout << " IA_energy : " << IA_energy << std::endl;

                double h_mu_mu = - IA_energy - ((Z_A - 0.5) * gamma_AA);
                double sum = 0;
                for (int B = 0; B < atoms_.size(); ++B) {
                    if (B != A) {
                        int Z_B = atoms_[B].getValElectrons();
                        sum += Z_B * gamma_matrix_(A, B);
                    }
                }
                H_core_matrix_(mu, nu) = h_mu_mu - sum;
            } else {
                // Off-diagonal elements (Equation 2.7)
            // std::cout << "(mu,nu) : (" << mu << "," << nu << ")" << std::endl;
                CNDO2Parameters params_mu = semiEmpiricalParams_.at(basisFunctions_[mu].element);
                CNDO2Parameters params_nu = semiEmpiricalParams_.at(basisFunctions_[nu].element);
            // std::cout << "Element_mu_beta = " << params_mu.beta << " & Element_nu_beta = " << params_nu.beta << std::endl;

                double beta_sum = (params_mu.beta + params_nu.beta);
                double s_uv = S_overlap_matrix_(mu, nu);
                H_core_matrix_(mu, nu) = 0.5 * beta_sum * s_uv;
            }
        }
    }
}

// Function to calculate the atomic electron populations from the basis function density matrix
void Molecule::calculate_pTot() {
    int num_atoms = atoms_.size();
    pTot_.set_size(num_atoms);
    pTot_.zeros(); // Set all elements to zero

    // Sum the densities from pTot_matrix_ for each atom
    for (int mu = 0; mu < N_; ++mu) {
        int atom_mu = basisFunctions_[mu].atomic_pos;  // Get the atom index for this basis function
        // Add the diagonal elements from the total density matrix for basis functions centered on atom_mu
        pTot_(atom_mu) += pA_matrix_(mu, mu) + pB_matrix_(mu, mu);
    }
}


// Calculates the total electronic energy of the molecule.
double Molecule::calc_totalE() {
    double conversion_factor = 27.211324570273; // a.u. to eV

    // Alpha and beta electron energies
    double alpha_electron_E_ = 0.5 * arma::accu(pA_matrix_ % (H_core_matrix_ + FA_matrix_));
    double beta_electron_E_ = 0.5 * arma::accu(pB_matrix_ % (H_core_matrix_ + FB_matrix_));

    // Nuclear repulsion energy
    for (int A = 0; A < atoms_.size(); ++A) {
        for (int B = A + 1; B < atoms_.size(); ++B) {
            int Z_A = atoms_[A].getValElectrons();
            int Z_B = atoms_[B].getValElectrons();
            arma::vec RA = arma::conv_to<arma::vec>::from(atoms_[A].getCenter());
            arma::vec RB = arma::conv_to<arma::vec>::from(atoms_[B].getCenter());
            double R_AB = arma::norm(RA - RB);
            nuclear_repulsion_E_ += (Z_A * Z_B) / R_AB;
        }
    }

    // Total energy
    return alpha_electron_E_ + beta_electron_E_ + (nuclear_repulsion_E_ * conversion_factor) ;
}


// // Constructs the orthogonalized X matrix for the molecule.
// void Molecule::build_X() {
//     // Initalize the eigval vector and eigvec matrix
//     arma::vec eigval;
//     arma::mat eigvec;
//     // Use Armadillo to solve system for S overlap matrix
//     arma::eig_sym(eigval, eigvec, S_overlap_matrix_);
//     // Use arma::pinv() to calculate the pseudo-inverse of sqrt(eigenvalues) matrix
//     arma::mat eig_inv_sqrt = arma::inv(arma::diagmat(arma::sqrt(eigval)));   
//     // Construct X
//     X_matrix_ = eigvec * eig_inv_sqrt * eigvec.t();
// }

// // Constructs the molecular orbital coefficient matrix (C_matrix_).
// void Molecule::build_C(){
//     // Form the Hamiltonian in the ortho basis
//     orthoH_matrix_ = X_matrix_.t() * H_matrix_ * X_matrix_;

//     // Diagonalize the orthogonalized Hamiltonian matrix using Armadillo
//     arma::eig_sym(epsilon_, V_matrix_, orthoH_matrix_);

//     // C = X * V 
//     C_matrix_ = X_matrix_ * V_matrix_;
// }

// // Constructs the molecular orbital overlap matrix (MO_S_matrix_).
// void Molecule::build_MO_S(){
//     // Form the MO_S with MO_coeffs and S
//     MO_S_matrix_ = C_matrix_.t() * S_overlap_matrix_ * C_matrix_;
// }

