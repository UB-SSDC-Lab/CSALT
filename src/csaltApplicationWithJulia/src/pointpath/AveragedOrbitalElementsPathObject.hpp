
#ifndef AVERAGED_ORBITAL_ELEMENTS_PATH_OBJECT_HPP
#define AVERAGED_ORBITAL_ELEMENTS_PATH_OBJECT_HPP

#include "UserPathFunction.hpp"
#include "julia_utils.hpp" 

class AveragedOrbitalElementsPathObject : public UserPathFunction 
{
public:
    // Default constructor
    AveragedOrbitalElementsPathObject(Integer n);

    // Destructor
    virtual ~AveragedOrbitalElementsPathObject();

    // Set Gravitational parameter
    void SetGravitationalParameter(Real muarg);

    // Set thrust parameters
    void SetThrustParameters(Real tMaxarg, Real Isparg, Real g0arg);

    // Evaluate functions
    void EvaluateFunctions();

    // Evaluate jacobians
    void EvaluateJacobians();

protected:

    // Dynamical parameters
    Real mu;
    Real tMax;
    Real Isp;
    Real g0;

    // State and derivative Julia vectors
    //JLvector x;
    //JLvector dx;

    // Gauss quadrature weights and notes
    Integer n;
    //JLvector taus;
    //JLvector ws;
};

#endif