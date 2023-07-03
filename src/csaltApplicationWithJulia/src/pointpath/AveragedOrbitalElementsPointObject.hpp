
#ifndef AVERAGED_ORBITAL_ELEMENTS_POINT_OBJECT_HPP
#define AVERAGED_ORBITAL_ELEMENTS_POINT_OBJECT_HPP

#include "UserPointFunction.hpp"

class AveragedOrbitalElementsPointObject : public UserPointFunction
{
public:
    // Default constructor
    AveragedOrbitalElementsPointObject();

    // Copy constructor
    AveragedOrbitalElementsPointObject(const AveragedOrbitalElementsPointObject& copy);

    // operator= 
    AveragedOrbitalElementsPointObject& operator=(const AveragedOrbitalElementsPointObject& copy);

    // Destructor
    virtual ~AveragedOrbitalElementsPointObject();

    // Set terminal state constraints
    void SetInitialStateConstraint(Rvector x0in);
    void SetFinalStateConstraint(Rvector xfin);

    // Evaluate functions
    void EvaluateFunctions();

    // Evaluate jacobians
    void EvaluateJacobians();

protected:
    // Initial time constraint
    Real t0_con;

    // Initial state constraint
    Rvector x0_con;

    // Final state constraint
    Rvector xf_con;
};

#endif