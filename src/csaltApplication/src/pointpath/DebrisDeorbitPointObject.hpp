
#ifndef DEBRIS_DEORBIT_POINT_OBJECT_HPP
#define DEBRIS_DEORBIT_POINT_OBJECT_HPP

#include "UserPointFunction.hpp"

class DebrisDeorbitPointObject : public UserPointFunction
{
public:
    // Default constructor
    DebrisDeorbitPointObject();

    // Copy constructor
    DebrisDeorbitPointObject(const DebrisDeorbitPointObject& copy);

    // operator=
    DebrisDeorbitPointObject& operator=(const DebrisDeorbitPointObject& copy);

    // Default destructor
    virtual ~DebrisDeorbitPointObject();

    // Set time constraints
    void SetTimeConstraints(Real t0arg, Real tfminarg, Real tfmaxarg);

    // Evaluate functions
    void EvaluateFunctions();

    // Evaluate jacobians
    void EvaluateJacobians();

protected:
    // Time constraints
    Real t0;
    Real tfmin;
    Real tfmax;
};

#endif