// Chem 279: Numerical Algorithms Applied to Computational Quantum Chemistry            
// Creator: Brandon Robello
// Date Created: 10/18/24

// ExtendedHuckel_test.cpp - this is an executable file that tests the 
// current state of the repo with input sample files

#include "AtomReader.h"
#include "Reader.h"
#include "Molecule.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>  
#include <filesystem>


// // Function to convert energy from eV to kJ/mol
// double convert_ev_to_kjmol(double energy_ev) {
//     // Conversion factor from eV to kJ/mol
//     const double conversion_factor = 96.485;

//     // Convert energy from eV to kJ/mol
//     return energy_ev * conversion_factor;
// }

int main(int argc, char* argv[]) {
    // Create an output file stream to redirect std::cout
    std::ofstream outFile("../Utils/Outputs/SCF_test_log.txt");  // File to store all cout output
    if (!outFile) {
        std::cerr << "Error: Unable to open output file for logging." << std::endl;
        return 1;
    }

    // Redirect std::cout to the output file
    std::streambuf* coutBuffer = std::cout.rdbuf();  // Save old buffer for restoration later
    std::cout.rdbuf(outFile.rdbuf());  // Redirect cout to outFile

    // Print out to terminal the start of the executable
    std::cout << "TESTING" << std::endl;

    // Set the directory path to the input files and display path that is being iterated in
    std::string inputDirPath = "../Utils/Inputs";
    std::cout << "Input Directory path is \""<< inputDirPath <<"\"" << std::endl;

    // Create a vector to store the individual input file paths
    std::vector<std::string> filePaths;

    // Iterate over the input directory and get the file paths and print the files
    std::cout << "Input Files:" << std::endl;
    for (const auto& entry : std::filesystem::directory_iterator(inputDirPath)) {
        if (entry.is_regular_file()) {  // Only process regular files
            std::filesystem::path filename = entry.path().filename(); // Convert path to filename to print
            filePaths.push_back(entry.path().string());  // Convert path to string and store it
            // Echo filenames found
            std::cout << filename << std::endl;
        }   
    }

    // Create reader object for expected format.
    AtomReader atomreader;

    // Loop over each file and pass it to the Reader function
    for (const auto& filePath : filePaths) {
        
        // Notify Processing has begun
        std::cout << "_______________________________________________________________" << std::endl;
        std::cout << "Processing file: " << filePath << std::endl;
        
        // Try-catch block in case the reader finds an atom not allowed
        try {
            // Call the AtomReader function for each file
            std::vector<Atom> atoms = atomreader.readFile(filePath);

            // Check if any Gaussians were inputted
            if (atoms.empty()) { 
                std::cout << "No atoms in system" << std::endl;
                continue;
            } 

            // Debugging IO
            // int iter = 1;
            // for (const auto& atom : atoms) {
            //     std::cout << "atom " << iter << ": " << atom << std::endl;
            //     ;
            //     iter++;
            // }

            Molecule molecule(atoms);
        
            // Print result
            std::cout.precision(5);
            std::cout << "_______________________________________________________________" << std::endl;
            std::cout << "Molecule: \n" << molecule << std::endl;

            molecule.calculate_gamma();

            molecule.printGamma();

            molecule.printS();

            molecule.build_core_hamiltonian();

            molecule.printH();

            molecule.build_Fmatrix_wSCF();

            double totalE = molecule.calc_totalE();

            std::cout << "Molecule specificied in " << filePath << " has an energy of " << totalE << " eV" << std::endl;

            // double total_energy = molecule.calc_totalE();

            // std::cout << "The molecule in file " << filePath << " has total energy: " << total_energy << std::endl;

            // if (molecule.getAtoms().size() == 2){
            //     // Energy molecules calculated
            //     double H_E = -13.6;
            //     // Difference of H2 from 2 Hs
            //     double bond_energy = total_energy - (2 * H_E);

            //     std::cout << "The bond energy of H2 is " << bond_energy << " eV." << std::endl;
            //     std::cout << "The difference between the total energy of H2 (" 
            //     << total_energy << " eV) and twice the energy of an H atom (" 
            //     << 2 * H_E << " eV)." << std::endl;
            // }
            // else if (molecule.getAtoms().size() == 6){
            //     // Energy molecules calculated
            //     double H2_E = -35.31;
            //     double C2H2_E = -177.17;

            //     // Difference of C2H4 and C2H2 + H2
            //     double reaction_energy_diff = total_energy - C2H2_E - H2_E;

            //     std::cout << "The energy difference for the reaction C2H2 + H2 -> C2H4 is " 
            //     << reaction_energy_diff << " eV or " << convert_ev_to_kjmol(reaction_energy_diff) << " Kj/mol."<< std::endl;
            //     std::cout << "The total energy of C2H4 (" << total_energy << " eV), minus the sum of C2H2 (" 
            //     << C2H2_E << " eV) and H2 (" << H2_E << " eV)." << std::endl;
            // }
            // else if (molecule.getAtoms().size() == 14){
            //     // Energy molecules calculated
            //     double H2_E = -35.31;
            //     double C6H6_E = -429.77;

            //     // Difference of C2H4 and C2H2 + H2
            //     double reaction_energy_diff = total_energy - C6H6_E - H2_E;

            //     std::cout << "The energy difference for the reaction C6H6 + H2 -> C6H8 is " 
            //     << reaction_energy_diff << " eV or " << convert_ev_to_kjmol(reaction_energy_diff) << " Kj/mol."<< std::endl;
            //     std::cout << "The total energy of C6H8 (" << total_energy << " eV), minus the sum of C6H6 (" 
            //     << C6H6_E << " eV) and H2 (" << H2_E << " eV)." << std::endl;
            // }
            

        } catch (const std::exception& e) {
            // Catch any exception and print the error message, but continue processing next file
            std::cerr << "Error reading file " << filePath << ": " << e.what() << std::endl;
            std::cerr << "Skipping file and continuing to the next...\n";
        }
    }

    // Notify that all files have been processed
    std::cout << "_______________________________________________________________" << std::endl;
    std::cout << "\nAll files processed.\n";

    // Restore std::cout to its original buffer (the terminal)
    std::cout.rdbuf(coutBuffer);

    // Close the output file
    outFile.close();

    return 0;
}
