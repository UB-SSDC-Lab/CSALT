
#ifndef DEBRIS_DEORBIT_PATH_OBJECT_HPP
#define DEBRIS_DEORBIT_PATH_OBJECT_HPP

#include "UserPathFunction.hpp"

class DebrisDeorbitPathObject : public UserPathFunction
{
public:
    // Default constructor
    DebrisDeorbitPathObject();

    // Copy constructor
    DebrisDeorbitPathObject(const DebrisDeorbitPathObject& copy);

    // operator=
    DebrisDeorbitPathObject& operator=(const DebrisDeorbitPathObject& copy);

    // Default destructor
    virtual ~DebrisDeorbitPathObject();

    // Set gravitational parameter
    void SetGravitationalParameter(Real muarg);

    // Set maximum thrust
    void SetMaximumThrust(Real tMaxarg);

    // Set mass parameters
    void SetMassParameters(Real m1arg, Real m2arg, Real mtarg);

    // Set tether parameters
    void SetTetherParameters(Real l0arg, Real karg, Real carg);

    // Evaluate functions
    void EvaluateFunctions();

    // Evaluate jacobians
    void EvaluateJacobians();

protected:
    // Compute dynamics partials
    Rmatrix ComputeStateDynamicsPartials();
    Rmatrix ComputeControlDynamicsPartials();
    Rmatrix ComputeTimeDynamicsPartials();

    // Compute path function partials
    Rmatrix ComputeStateConstraintPartials();
    Rmatrix ComputeControlConstraintPartials();
    Rmatrix ComputeTimeConstraintPartials();

    // Gravitational parameter
    Real mu;

    // Spacecraft's maximum thrust
    Real tMax;

    // Orbital inclination (not required for two-dim model)
    Real inc;

    // Mass parameters
    Real m1;
    Real m2;
    Real mt;
    Real m;
    Real ms;

    // Tether parameters 
    Real l0;
    Real k;
    Real c;
};

#endif