// Chem279: Numerical Algorithms Applied to Computational Quantum Chemistry  
// Creator: Brandon Robello
// Date Created: 10/2/24
//
// Reader.h contains the necessary include files and 
//        C++ Class declarations for Reader functions

#ifndef READER_H
#define READER_H

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>  


/**
 * @class Reader
 * @brief Abstract base class for reading objects from files.
 */
template<typename T>
class Reader {
public:
    /**
     * @brief Reads objects from the specified file.
     * 
     * @param filename The name of the file to read from.
     * @return A vector of Gaussian objects.
     */
    // Base class function to read a file and delegate line parsing to derived classes
    std::vector<T> readFile(const std::string& filename) const {
        std::vector<T> objects;
        std::ifstream file(filename);  // Open the file

        if (!file) {  // Check if the file was opened successfully
            std::cerr << "Error: Unable to open the file " << filename << std::endl;
            return objects;  // Return empty vector if file couldn't be opened
        }

        // Read the first line to get the expected number of lines
        int num_lines = 0;
        if (!(file >> num_lines)) {
            std::cerr << "Error: Could not read the number of lines from the first line" << std::endl;
            return objects;
        }

        std::string line;  // String to hold each line of the file
        std::getline(file, line); // Move to the next line after reading num_lines

        // Read exactly `num_lines` lines
        for (int i = 0; i < num_lines; ++i) {
            if (!std::getline(file, line)) {
                std::cerr << "Error: Reached end of file before reading the expected number of lines" << std::endl;
                break;
            }

            try {
                T obj = parseLine(line);  // Parse each line
                objects.emplace_back(obj);
            } catch (const std::invalid_argument& e) {
                std::cerr << "Error processing line " << i + 1 << ": " << line << "\n" << e.what() << std::endl;
            }
        }

        // Check if the actual number of lines read matches the expected number
        if (objects.size() != num_lines) {
            std::cerr << "Warning: Expected " << num_lines << " atoms but only read " << objects.size() << std::endl;
        }

        file.close();  // Close the file
        return objects;  // Return the vector of objects
    }

protected:
    /**
     * @brief Parses a line and creates a  object.
     *        This function must be implemented by derived classes.
     * 
     * @param line The line to be parsed.
     * @return A object created from the line.
     */
    virtual T parseLine(const std::string& line) const = 0;

    virtual ~Reader() = default; // Virtual destructor for proper cleanup of derived classes
};


#endif
