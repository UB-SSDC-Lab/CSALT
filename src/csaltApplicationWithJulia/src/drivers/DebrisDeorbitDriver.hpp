
#ifndef DEBRIS_DEORBIT_DRIVER_HPP
#define DEBRIS_DEORBIT_DRIVER_HPP

#include "CsaltDriver.hpp"

class DebrisDeorbitDriver : public CsaltDriver
{
public:
    DebrisDeorbitDriver();
    DebrisDeorbitDriver(const std::string &testName, const Integer snoptConsoleOutputLevel);
    virtual ~DebrisDeorbitDriver();

    // Process trajectory solution
    virtual void ProcessSolution(Trajectory* traj);

protected:
    virtual void SetParameters();
    virtual void ProcessInputFile();
    virtual void SetPointPathAndProperties();
    virtual void SetOptimalControlFunctions(std::vector<OptimalControlFunction*>& funcList);
    virtual void SetupPhases();

    // Meta information
    std::string mode;
    Integer maxMeshRefinements;
    Integer maxMajorIterations;
    Integer maxTotalIterations;
    Integer maxPolyDegIncrease;
    Real majorOptimalityTol;
    Real feasibilityTol;
    Real meshRelTol;

    // Guess information
    std::string ochFile;
    Real guess_t0;
    Real guess_tf;

    // Parameters
    Real mu;
    Real tMax;
    Real m1;
    Real m2;
    Real mt;
    Real l0;
    Real k;
    Real c;

    // Problem units (not currently used)
    Real TU;
    Real LU;

    // Bounds
    Real time_LB;
    Real time_UB;
    Rvector state_LB;
    Rvector state_UB;
    Rvector control_LB;
    Rvector control_UB;

    // Mesh
    Rvector meshFractions;
    IntegerArray meshNumPoints;

    // Initial state constraints
    Real t0;
    Real a0;
    bool a0_con_on;
    Real e0;
    bool e0_con_on;
    Real aop0;
    bool aop0_con_on;
    Real ta0;
    bool ta0_con_on;
    Real alpha0;
    bool alpha0_con_on;
    Real d_alpha0;
    bool d_alpha0_con_on;
    Real L0;
    bool L0_con_on;
    Real d_L0;
    bool d_L0_con_on;

    // Final state constraints
    Real rpf;
    bool rpf_con_on;
    Real tfmin;
    Real tfmax;

    // Define an enum for reading input file
    enum InputFileLocation
    {
        NS,
        META,
        GUESS, 
        PARAMETERS,
        BOUNDS,
        MESH,
        CONSTRAINTS
    };
};

#endif