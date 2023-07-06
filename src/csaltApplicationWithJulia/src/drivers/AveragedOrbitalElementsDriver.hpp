
#ifndef AVERAGED_ORBITAL_ELEMENTS_DRIVER_HPP
#define AVERAGED_ORBITAL_ELEMENTS_DRIVER_HPP

#include "CsaltDriver.hpp"
#include "julia.h"

class AveragedOrbitalElementsDriver : public CsaltDriver 
{
public:
    AveragedOrbitalElementsDriver();
    AveragedOrbitalElementsDriver(const std::string &testName, const Integer snoptConsoleOutputLevel);
    virtual ~AveragedOrbitalElementsDriver();

    // Process trajectory solution
    virtual void ProcessSolution(Trajectory* traj);

protected:
    virtual void SetParameters();
    virtual void SetPointPathAndProperties();
    virtual void SetupPhases();

    // Units
    Real TU;
    Real LU;

    // Gravitational parameter
    Real mu;

    // Problem parameters
    Real m0; 
    Real P;
    Real eta;
    Real Isp;  
    Real g0;   
    Real tMax;

    // Number of points used in Gauss quadrature
    Integer n;

    // Initial and final equality constrained states
    Rvector x0_con; // Should have length 6 (includes mass)
    Rvector xf_con; // Should have length 5 (does not include mass)
};

#endif