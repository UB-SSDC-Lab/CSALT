
#ifndef AVERAGED_ORBITAL_ELEMENTS_PATH_OBJECT_HPP
#define AVERAGED_ORBITAL_ELEMENTS_PATH_OBJECT_HPP

#include "UserPathFunction.hpp"
#include "csalt.hpp"
#include "julia_utils.hpp" 

using namespace jluna;

class AveragedOrbitalElementsPathObject : public UserPathFunction 
{
public:
    // Default constructor
    AveragedOrbitalElementsPathObject();

    // Destructor
    virtual ~AveragedOrbitalElementsPathObject();

    // Initialize path function
   virtual void Initialize(FunctionInputData *pd, PathFunctionContainer *pfc);

    // Set Gravitational parameter
    void SetGravitationalParameter(Real muarg);

    // Set thrust parameters
    void SetThrustParameters(Real tMaxarg, Real Isparg, Real g0arg);

    // Set gauss quadrature
    void SetGaussQuadrature(Integer narg);

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

    // Gauss quadrature weights and notes
    Integer n;
    Vector<Float64> taus;
    Vector<Float64> ws;
};

#endif