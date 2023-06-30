
#ifndef JULIA_UTILS_HPP
#define JULIA_UTILS_HPP


// A lot of code was taken from https://blog.esciencecenter.nl/10-examples-of-embedding-julia-in-c-c-66282477e62c#cc9a 

#include <string>
#include "julia.h"
//#include "csalt.hpp"

// Define constant string to julia base src directory
const std::string JL_SRC((char*) JLSRC);

// Starting and stoping julia (must call at start and end of program!)
void setup_julia();
void stop_julia(int n);

// Exception handling
void handle_julia_exception(void);
jl_value_t* handle_eval_string(const char* code);

// Vectors and matricies
class JLvector 
{
public:
    // Default constructor
    JLvector(size_t n);

    // Operators
    double&         operator()(size_t index) {return vec[index];}
    const double&   operator()(size_t index) const {return vec[index];}
    double&         operator[](size_t index) {return vec[index];}
    const double&   operator[](size_t index) const {return vec[index];}

    // Getters for raw data
    jl_array_t* GetJuliaVector() {return vecptr;}
    double* GetDataVector() {return vec;}

    // Make zero vector
    void MakeZeroVector();

protected:
    // Pointer to Julia vector
    jl_array_t* vecptr;

    // Pointer to vector
    double* vec;

    // Length of vector
    size_t n;
};

#endif