// Chem279: Numerical Algorithms Applied to Computational Quantum Chemistry  
// Creator: Brandon Robello
// Date Created: 9/18/24
// Adapted: 10/18/24
//          Added N_, Basis Functions
//          electrons_, valence electrons
//
// Atom.h contains the necessary include files and 
//        C++ Class declarations for Atom class

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
#include <set>
#include <map>

#ifndef ATOM_H
#define ATOM_H

/**
 * @class Atom
 * @brief Represents an atom with its atomic number and coordinates in 3D space.
 */
class Atom {
public:
    /**
     * @brief Constructs an Atom object with the specified atomic number and coordinates.
     * @param atomicNumber The atomic number of the element (e.g., 79 for Au).
     * @param x The x-coordinate of the atom in space.
     * @param y The y-coordinate of the atom in space.
     * @param z The z-coordinate of the atom in space.
     */
    Atom(int atomicNumber, double x, double y, double z);

    /**
     * @brief Overloads the << operator for printing the Atom object in a readable format.
     * @param os The output stream.
     * @param atom The Atom object to be output.
     * @return std::ostream& Reference to the output stream.
     */
    friend std::ostream& operator<<(std::ostream& os, const Atom& atom);

    /**
     * @brief Calculates the Euclidean distance between this atom and another atom.
     * @param other The other Atom object to measure the distance to.
     * @return double The Euclidean distance between the two atoms.
     */
    double distanceTo(const Atom& other) const;

    /**
     * @brief Calculates the difference in the x-coordinate between this atom and another atom.
     * @param other The other Atom object.
     * @return double The difference in the x-coordinate.
     */
    double get_dx(const Atom& other) const;

    /**
     * @brief Calculates the difference in the y-coordinate between this atom and another atom.
     * @param other The other Atom object.
     * @return double The difference in the y-coordinate.
     */
    double get_dy(const Atom& other) const;

    /**
     * @brief Calculates the difference in the z-coordinate between this atom and another atom.
     * @param other The other Atom object.
     * @return double The difference in the z-coordinate.
     */
    double get_dz(const Atom& other) const;

    /**
     * @brief Retrieves the N, number of basis functions for this Atom.
     * @return int The N of the element.
     */
    int getN() const;

    /**
     * @brief Retrieves the number of valence electrons for this Atom.
     * @return int The electrons of the element.
     */
    int getValElectrons() const;   

    /**
     * @brief Retrieves the atomic number of this Atom.
     * @return int The atomic number of the element.
     */
    int getAtomicNumber() const;

    /**
     * @brief Gets the center (x,y,z) position of the Atom
     */
    std::vector<double> getCenter() const;

    /**
     * @brief Retrieves the x-coordinate of this Atom.
     * @return double The x-coordinate.
     */
    double getX() const;

    /**
     * @brief Retrieves the y-coordinate of this Atom.
     * @return double The y-coordinate.
     */
    double getY() const;

    /**
     * @brief Retrieves the z-coordinate of this Atom.
     * @return double The z-coordinate.
     */
    double getZ() const;

    /**
     * @brief Sets a new x-coordinate for the Atom.
     * @param newX The new x-coordinate.
     */
    void setX(double newX);

    /**
     * @brief Sets a new y-coordinate for the Atom.
     * @param newY The new y-coordinate.
     */
    void setY(double newY);

    /**
     * @brief Sets a new z-coordinate for the Atom.
     * @param newZ The new z-coordinate.
     */
    void setZ(double newZ);

private:
    int ValElectrons_; ///< Number of valence electrons in atom.
    int N_;            ///< Number of functions for atom.
    int atomicNumber_; ///< The atomic number of the atom.
    double x_;         ///< The x-coordinate of the atom in space.
    double y_;         ///< The y-coordinate of the atom in space.
    double z_;         ///< The z-coordinate of the atom in space.
};

#endif
