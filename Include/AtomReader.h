// Chem279: Numerical Algorithms Applied to Computational Quantum Chemistry  
// Creator: Brandon Robello
// Date Created: 10/18/24
//
// AtomReader.h contains the necessary include files and 
//        C++ Class declarations for AtomReader functions

#ifndef ATOMREADER_H
#define ATOMREADER_H

#include "Reader.h"
#include "Atom.h"
#include <string>
#include <vector>

/**
 * @class AtomReader
 * @brief Derived class that reads Atom objects in the standard format (E, X, Y, Z).
 */
class AtomReader : public Reader<Atom> {
protected:
    /**
     * @brief Parses a line into Atom objects in the standard format (E, X, Y, Z).
     * 
     * @param line The line to be parsed.
     * @return A Atom object created from the line.
     */
    Atom parseLine(const std::string& line) const override;
};

#endif 