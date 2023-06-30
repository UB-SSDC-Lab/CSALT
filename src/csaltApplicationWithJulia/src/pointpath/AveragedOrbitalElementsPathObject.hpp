
#ifndef AVERAGED_ORBITAL_ELEMENTS_PATH_OBJECT_HPP
#define AVERAGED_ORBITAL_ELEMENTS_PATH_OBJECT_HPP

#include "UserPathFunction.hpp"
#include "julia_utils.hpp" 

class AveragedOrbitalElementsPathObject : public UserPathFunction 
{
    public:
        // Default constructor
        AveragedOrbitalElementsPathObject(Integer n);

        // Copy constructor
        AveragedOrbitalElementsPathObject(const AveragedOrbitalElementsPathObject& copy);

        // Assignment operator
        AveragedOrbitalElementsPathObject& operator=(const AveragedOrbitalElementsPathObject& copy);

        // Destructor
        virtual ~AveragedOrbitalElementsPathObject();

        // Set Gravitational parameter
        void SetGravitationalParameter(Real muarg);

        // Set thrust parameters
        void SetThrustParameters(Real tMax, Real Isp, Real g0);

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
    JLvector taus;
    JLvector ws;
};

#endif