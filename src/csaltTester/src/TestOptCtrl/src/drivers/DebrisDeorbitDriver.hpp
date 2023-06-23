
#ifndef DEBRIS_DEORBIT_DRIVER_HPP
#define DEBRIS_DEORBIT_DRIVER_HPP

#include "CsaltTestDriver.hpp"

class DebrisDeorbitDriver : public CsaltTestDriver
{
public:
    DebrisDeorbitDriver();
    DebrisDeorbitDriver(const std::string &testName, const Integer snoptConsoleOutputLevel);
    virtual ~DebrisDeorbitDriver();
protected:
    virtual void SetParameters();
    virtual void SetPointPathAndProperties();
    virtual void SetOptimalControlFunctions(std::vector<OptimalControlFunction*>& funcList);
    virtual void SetupPhases();

    // Problem units
    Real TU;
    Real LU;

    // Parameter bounds
    Real time_LB;
    Real time_UB;
    Rvector state_LB;
    Rvector state_UB;
    Rvector control_LB;
    Rvector control_UB;

    // Initial guess parameters
    Real guess_t0;
    Real guess_tf;
    std::string ochFile;

    // Path function constants
    Real mu;
    Real tMax;
    Real m1;
    Real m2;
    Real mt;
    Real l0;
    Real k;
    Real c;

    // Point function constants
    Real t0;
    Real tfmin;
    Real tfmax;
    Real a0;
    Real e0;
    Real aop0;
    Real ta0;
    Real alpha0;
    Real rpf;
};

#endif