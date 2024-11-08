# Using CMake version 3.14 for NumPy support
cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

# Specify project name & programming languages
project(HW3_ExtendedHuckel LANGUAGES CXX)

# Require a modern C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add Armadillo support
find_package(Armadillo REQUIRED)

# Add your existing include directories
include_directories(${PROJECT_SOURCE_DIR}/Include)
include_directories(${ARMADILLO_INCLUDE_DIRS})

# Add libraries built from Source/*.cpp
add_library(Atom ./Source/Atom.cpp)
add_library(AtomReader ./Source/AtomReader.cpp)
add_library(Molecule ./Source/Molecule.cpp)

# Add executables for tests
add_executable(SCF_test ./Test/SCF_test.cpp)

# Link libraries to tests and include necessary directories
target_include_directories(Atom PRIVATE ./Include)
target_include_directories(AtomReader PRIVATE ./Include)
target_include_directories(Molecule PRIVATE ./Include)

# Ensure Armadillo is linked with these targets
target_link_libraries(Molecule PRIVATE ${ARMADILLO_LIBRARIES})
target_link_libraries(AtomReader PRIVATE Atom )


# Link executables to their respective libraries and Armadillo
target_link_libraries(SCF_test PRIVATE AtomReader Molecule)
