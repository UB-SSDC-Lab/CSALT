
#include "DebrisDeorbitDriver.hpp"
#include "DebrisDeorbitPathObject.hpp"
#include "DebrisDeorbitPointObject.hpp"
#include "DebrisDeorbitInitialStateConstraint.hpp"
#include "DebrisDeorbitFinalStateConstraint.hpp"

DebrisDeorbitDriver::DebrisDeorbitDriver() : CsaltDriver("DebrisDeorbit") 
{
    SetParameters();
}

DebrisDeorbitDriver::DebrisDeorbitDriver(const std::string &testName, const Integer snoptConsoleOutputLevel) :
    CsaltDriver(testName, snoptConsoleOutputLevel) 
{
    SetParameters();
}

DebrisDeorbitDriver::~DebrisDeorbitDriver() {}

void DebrisDeorbitDriver::SetParameters()
{
    // Set problem units
    TU      = 1.0;
    LU      = 1.0;

    // Set parameter bounds
    time_LB = 0.0;
    time_UB = 10.0 * 3600.0;

    state_LB    = Rvector(8, 6000.0, 
                             1.0e-8, 
                             -GmatMathConstants::PI, 
                             0.0, 
                             -2.0*GmatMathConstants::PI, 
                             -0.5 * GmatMathConstants::PI / 180.0, 
                             30.0001e-3, 
                             -1.0e-3);
    state_UB    = Rvector(8, 10000.0,
                             0.1,
                             GmatMathConstants::PI,
                             100.0 * GmatMathConstants::PI,
                             2.0*GmatMathConstants::PI,
                             0.5 * GmatMathConstants::PI / 180.0,
                             35.0e-3,
                             1.0e-3);

    control_LB  = Rvector(3, 0.0, -1.0, -1.0);
    control_UB  = Rvector(3, 1.0,  1.0,  1.0);

    // Set guess parameters
    guess_t0    = 0.0;
    guess_tf    = 5.305479356551708 * 3600.0;
    ochFile     = "/Users/granthec/Documents/Projects/SCvxDebrisDeorbit/data/och_guess_files/guess.och";

    // Set gravitational parameter
    mu      = 3.986e5;

    // Set maximum thrust
    tMax    = 50.0e-3;

    // Set mass parameters
    m1      = 1600.0;
    m2      = 3000.0;
    mt      = 0.5;

    // Set tether parameters
    l0      = 30.0e-3;
    k       = 1573.0;
    c       = 16.0;

    // Set initial and final time constraints
    t0      = 0.0;
    tfmin   = 0.5  * 3600.0;
    tfmax   = 10.0 * 3600.0;

    // Set initial state constraints
    a0      = 6871.0;
    e0      = 0.001;
    aop0    = GmatMathConstants::PI / 2.0;
    ta0     = 0.0;
    alpha0  = -GmatMathConstants::PI / 2.0 + 0.01;

    // Set final periapsis constraint
    rpf     = 6478.0;
}

void DebrisDeorbitDriver::SetPointPathAndProperties()
{
    // Instantiate point and path objects
    pathObject = new DebrisDeorbitPathObject();
    pointObject = new DebrisDeorbitPointObject();

    // Set path object properties
    dynamic_cast<DebrisDeorbitPathObject*>(pathObject)->SetGravitationalParameter(mu);
    dynamic_cast<DebrisDeorbitPathObject*>(pathObject)->SetMaximumThrust(tMax);
    dynamic_cast<DebrisDeorbitPathObject*>(pathObject)->SetMassParameters(m1,m2,mt);
    dynamic_cast<DebrisDeorbitPathObject*>(pathObject)->SetTetherParameters(l0,k,c);

    // Set point object properties
    dynamic_cast<DebrisDeorbitPointObject*>(pointObject)->SetTimeConstraints(t0,tfmin,tfmax);

    // Instantiate vector of optimal control functions
    std::vector<OptimalControlFunction*> funcList;
    SetOptimalControlFunctions(funcList);
    pointObject->AddFunctions(funcList);

    // Set properties
    maxMeshRefinementCount      = 20;
    majorIterationsLimits[0]    = 500;
    totalIterationsLimits[0]    = 200000;
    optimizationMode            = StringArray(1, "Minimize");

    majorOptimalityTolerances.SetSize(1);
    feasibilityTolerances.SetSize(1);
    majorOptimalityTolerances(0)    = 1e-5;
    feasibilityTolerances(0)        = 1e-5;
}

void DebrisDeorbitDriver::SetOptimalControlFunctions(std::vector<OptimalControlFunction*>& funcList)
{
    bool initStateCon   = true;
    bool finStateCon    = true;

    // ===== Initial state constraint
    if (initStateCon) {
        // Set state constraints
        Rvector initialState(5, a0, e0, aop0, ta0, alpha0);

        // Set state constraint idxs
        IntegerArray initialStateIdxs;
        for (Integer i = 0; i < 5; i++)
            initialStateIdxs.push_back(i);

        // Set upper and lower bounds
        Rvector functionLB(5);
        Rvector functionUB(5);
        functionLB.MakeZeroVector();
        functionUB.MakeZeroVector();

        // Define phase dependencies
        IntegerArray phaseDepends(1);
        phaseDepends.at(0) = 0;

        // Define point dependencies
        IntegerArray pointDepends(1);
        pointDepends.at(0) = 0;

        // Instantiate constraint
        DebrisDeorbitInitialStateConstraint* initCon = new DebrisDeorbitInitialStateConstraint(
            "initial_state_constraint");

        // Initialize constraint
        initCon->Initialize();
        initCon->SetNumFunctions(5);
        initCon->SetPhaseList(phaseList);
        initCon->SetPhaseDependencies(phaseDepends);
        initCon->SetPointDependencies(pointDepends);
        initCon->SetLowerBounds(functionLB);
        initCon->SetUpperBounds(functionUB);
        initCon->SetInitialStateConstraint(initialState);
        initCon->SetInitialStateConstraintIdxs(initialStateIdxs);
        initCon->SetAnalyticJacobianOn(true);
        funcList.push_back(initCon);
    }

    // Final state constraint
    if (finStateCon) {
        // Set upper and lower bounds
        Rvector functionLB(1);
        Rvector functionUB(1);
        functionLB.MakeZeroVector();
        functionUB.MakeZeroVector();

        // Define phase dependencies
        IntegerArray phaseDepends(1);
        phaseDepends.at(0) = 0;

        // Define point dependencies
        IntegerArray pointDepends(1);
        pointDepends.at(0) = 1;

        // Instantiate constraint
        DebrisDeorbitFinalStateConstraint* finCon = new DebrisDeorbitFinalStateConstraint(
            "final_state_constraint");

        // Initialize constraint
        finCon->Initialize();
        finCon->SetNumFunctions(1);
        finCon->SetPhaseList(phaseList);
        finCon->SetPhaseDependencies(phaseDepends);
        finCon->SetPointDependencies(pointDepends);
        finCon->SetLowerBounds(functionLB);
        finCon->SetUpperBounds(functionUB);
        finCon->SetFinalStateConstraint(rpf);
        finCon->SetAnalyticJacobianOn(true);
        funcList.push_back(finCon);
    }
}

void DebrisDeorbitDriver::SetupPhases()
{
    // Instantiate LGR phase
    RadauPhase *phase1              = new RadauPhase();

    // Set guess mode
    std::string initialGuessMode    = "OCHFile";

    // Set mesh properties
    Integer n = 10;
    Integer m = 12;
    Real step = 2.0 / (n - 1);
    Rvector meshIntervalFractions(n);
    IntegerArray meshIntervalNumPoints;
    for (Integer i = 0; i < n - 1; i++) {
        meshIntervalFractions(i) = -1.0 + step*i;
        meshIntervalNumPoints.push_back(m);
    }
    meshIntervalFractions(n - 1) = 1.0;

    // Set phase properties 
    phase1->SetRelativeErrorTol(5e-6);
    phase1->SetInitialGuessMode(initialGuessMode);
    phase1->SetGuessFileName(ochFile);
    phase1->SetNumStateVars(8);
    phase1->SetNumControlVars(3);
    phase1->SetMeshIntervalFractions(meshIntervalFractions);
    phase1->SetMeshIntervalNumPoints(meshIntervalNumPoints);
    phase1->SetStateLowerBound(state_LB);
    phase1->SetStateUpperBound(state_UB);
    phase1->SetControlLowerBound(control_LB);
    phase1->SetControlUpperBound(control_UB);
    phase1->SetTimeLowerBound(time_LB);
    phase1->SetTimeUpperBound(time_UB);
    phase1->SetTimeInitialGuess(guess_t0);
    phase1->SetTimeFinalGuess(guess_tf);

    // Push phase to phase list
    phaseList.push_back(phase1);
}