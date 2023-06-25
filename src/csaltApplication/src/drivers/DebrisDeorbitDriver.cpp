
#include "DebrisDeorbitDriver.hpp"
#include "DebrisDeorbitPathObject.hpp"
#include "DebrisDeorbitPointObject.hpp"
#include "DebrisDeorbitInitialStateConstraint.hpp"
#include "DebrisDeorbitFinalStateConstraint.hpp"

#include <fstream>
#include <sstream>
#include <string>

DebrisDeorbitDriver::DebrisDeorbitDriver() : CsaltDriver("DebrisDeorbit") {}

DebrisDeorbitDriver::DebrisDeorbitDriver(const std::string &testName, const Integer snoptConsoleOutputLevel) :
    CsaltDriver(testName, snoptConsoleOutputLevel) {}

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
                             -5.0 * GmatMathConstants::PI / 180.0, 
                             30.0001e-3, 
                             -1.0e-3);
    state_UB    = Rvector(8, 10000.0,
                             0.1,
                             GmatMathConstants::PI,
                             100.0 * GmatMathConstants::PI,
                             2.0*GmatMathConstants::PI,
                             5.0 * GmatMathConstants::PI / 180.0,
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
    // Process the input file 
    ProcessInputFile();

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
    maxMeshRefinementCount      = maxMeshRefinements;
    majorIterationsLimits[0]    = maxMajorIterations;
    totalIterationsLimits[0]    = maxTotalIterations;
    if (mode == "MINIMIZE")
        optimizationMode.push_back("Minimize");
    else if (mode == "FP")
        optimizationMode.push_back("Feasible point");
    else {
        std::cout << "Unsupported option for MODE, setting mode to FP";
        optimizationMode.push_back("Feasible point");
    }

    majorOptimalityTolerances.SetSize(1);
    feasibilityTolerances.SetSize(1);
    majorOptimalityTolerances(0)    = majorOptimalityTol;
    feasibilityTolerances(0)        = feasibilityTol;
}

void DebrisDeorbitDriver::SetOptimalControlFunctions(std::vector<OptimalControlFunction*>& funcList)
{
    bool initStateCon   = true;
    bool finStateCon    = true;

    // ===== Initial state constraint
    if (initStateCon) {
        // Allocate vector we can push constraints and their indecies to
        std::vector<Real> initialStateVec;
        IntegerArray initialStateIdxs;

        // Determing the number of initial state constraints
        if (a0_con_on) {
            initialStateVec.push_back(a0);
            initialStateIdxs.push_back(0);
        }
        if (e0_con_on) {
            initialStateVec.push_back(e0);
            initialStateIdxs.push_back(1);
        }
        if (aop0_con_on) {
            initialStateVec.push_back(aop0);
            initialStateIdxs.push_back(2);
        }
        if (ta0_con_on) {
            initialStateVec.push_back(ta0);
            initialStateIdxs.push_back(3);
        }
        if (alpha0_con_on) {
            initialStateVec.push_back(alpha0);
            initialStateIdxs.push_back(4);
        }
        if (d_alpha0_con_on) {
            initialStateVec.push_back(d_alpha0);
            initialStateIdxs.push_back(5);
        }
        if (L0_con_on) {
            initialStateVec.push_back(L0);
            initialStateIdxs.push_back(6);
        }
        if (d_L0_con_on) {
            initialStateVec.push_back(d_L0);
            initialStateIdxs.push_back(7);
        }

        // Set state constraints
        Rvector initialState(initialStateVec.size());
        for (Integer i = 0; i < initialStateVec.size(); i++)
            initialState(i) = initialStateVec.at(i);

        // Set upper and lower bounds
        Rvector functionLB(initialStateVec.size());
        Rvector functionUB(initialStateVec.size());
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
        initCon->SetNumFunctions(initialStateVec.size());
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
    if (finStateCon && rpf_con_on) {
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
    Rvector meshIntervalFractions;
    IntegerArray meshIntervalNumPoints;
    if (meshFractions.GetSize() > 1) {
        // If mesh properties directly set in input file, just set
        meshIntervalFractions = meshFractions;
        meshIntervalNumPoints = meshNumPoints;
    }
    else {
        // Only number of mesh fractions and points specified, create mesh
        Integer n = Integer (meshFractions(0));
        Integer m = meshNumPoints[0];
        Real step = 2.0 / (meshFractions(0) - 1.0);
        meshIntervalFractions.SetSize(n);
        for (Integer i = 0; i < n - 1; i++) {
            meshIntervalFractions(i) = -1.0 + step*i;
            meshIntervalNumPoints.push_back(m);
        }
        meshIntervalFractions(n - 1) = 1.0;
    }

    // Set phase properties 
    phase1->SetRelativeErrorTol(meshRelTol);
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

void DebrisDeorbitDriver::ProcessInputFile()
{
    // Instantiate file stream
    std::ifstream infile(inputFile);

    // Loop through file line by line
    InputFileLocation loc = NS;
    std::string line, key, val;
    while (std::getline(infile, line))
    {
        if (loc == NS) { // Not in section currently, see if we're entering a section
            if (line.find("META_START") != std::string::npos)
                loc = META;
            else if (line.find("GUESS_START") != std::string::npos)
                loc = GUESS;
            else if (line.find("PARAMETERS_START") != std::string::npos)
                loc = PARAMETERS;
            else if (line.find("BOUNDS_START") != std::string::npos) 
                loc = BOUNDS;
            else if (line.find("MESH_START") != std::string::npos)
                loc = MESH;
            else if (line.find("CONSTRAINTS_START") != std::string::npos)
                loc = CONSTRAINTS;
        }
        else {
            if (loc == META) {
                if (line.find("META_STOP") != std::string::npos)
                    loc = NS;
                else { // Process meta data
                    // Grab key and value
                    GetKeyValuePair(line, key, val);

                    // Set value for key
                    if (key.find("MODE") != std::string::npos)
                        mode = val;                    
                    if (key.find("SOL_FILE") != std::string::npos) {
                        controlHistoryFile = val;
                        controlHistoryFile.erase(std::remove(controlHistoryFile.begin(),
                                                             controlHistoryFile.end(),
                                                             '\"'),
                                                 controlHistoryFile.end());
                    }
                    else if (key.find("MAX_MESH_REFINEMENTS") != std::string::npos)
                        maxMeshRefinements = std::stoi(val); 
                    else if (key.find("MAX_MAJOR_ITERATIONS") != std::string::npos)
                        maxMajorIterations = std::stoi(val);
                    else if (key.find("MAX_TOTAL_ITERATIONS") != std::string::npos)
                        maxTotalIterations = std::stoi(val);
                    else if (key.find("MAJOR_OPTIMALITY_TOL") != std::string::npos)
                        majorOptimalityTol = std::stod(val);
                    else if (key.find("FEASIBILITY_TOL") != std::string::npos)
                        feasibilityTol = std::stod(val);
                    else if (key.find("MESH_REL_TOL") != std::string::npos) 
                        meshRelTol = std::stod(val);
                }
            }
            else if (loc == GUESS) {
                if (line.find("GUESS_STOP") != std::string::npos)
                    loc = NS;
                else { // Process guess data
                    // Grab key and value
                    GetKeyValuePair(line, key, val);

                    // Set value for key
                    if (key.find("GUESS_FILE") != std::string::npos) {
                        ochFile = val;
                        ochFile.erase(std::remove(ochFile.begin(),ochFile.end(),'\"'),ochFile.end());
                    }
                    else if (key.find("GUESS_T0") != std::string::npos)
                        guess_t0 = std::stod(val);
                    else if (key.find("GUESS_TF") != std::string::npos)
                        guess_tf = std::stod(val);
                }
            }
            else if (loc == PARAMETERS) {
                if (line.find("PARAMETERS_STOP") != std::string::npos)
                    loc = NS;
                else { // Process parameter data
                    // Grab key and value
                    GetKeyValuePair(line, key, val);

                    // Set value for key
                    if (key.find("MU") != std::string::npos)
                        mu = std::stod(val);
                    else if (key.find("T_MAX") != std::string::npos)
                        tMax = std::stod(val);
                    else if (key.find("M1") != std::string::npos)
                        m1 = std::stod(val);
                    else if (key.find("M2") != std::string::npos)
                        m2 = std::stod(val);
                    else if (key.find("MT") != std::string::npos)
                        mt = std::stod(val);
                    else if (key.find("L0") != std::string::npos)
                        l0 = std::stod(val);
                    else if (key.find("K") != std::string::npos)
                        k = std::stod(val);
                    else if (key.find("C") != std::string::npos)
                        c = std::stod(val);
                }
            }
            else if (loc == BOUNDS) {
                if (line.find("BOUNDS_STOP") != std::string::npos)
                    loc = NS;
                else { // Process bounds data
                    // Grab key and value
                    GetKeyValuePair(line, key, val);

                    // Set value for key
                    if (key.find("TIME_LB") != std::string::npos)
                        time_LB = std::stod(val);
                    else if (key.find("TIME_UB") != std::string::npos)
                        time_UB = std::stod(val);
                    else if (key.find("STATE_LB") != std::string::npos)
                        state_LB = StringToRvector(val);
                    else if (key.find("STATE_UB") != std::string::npos)
                        state_UB = StringToRvector(val);
                    else if (key.find("CONTROL_LB") != std::string::npos)
                        control_LB = StringToRvector(val);
                    else if (key.find("CONTROL_UB") != std::string::npos)
                        control_UB = StringToRvector(val);
                }
            }
            else if (loc == MESH) {
                if (line.find("MESH_STOP") != std::string::npos)
                    loc = NS;
                else { // Process mesh data
                    // Grab key and value
                    GetKeyValuePair(line, key, val);

                    // Set value for key
                    if (key.find("MESH_FRACS") != std::string::npos)
                        meshFractions = StringToRvector(val);
                    else if (key.find("MESH_NPOINT") != std::string::npos)
                        meshNumPoints = StringToIntegerArray(val);
                }
            }
            else if (loc == CONSTRAINTS) {
                if (line.find("CONSTRAINTS_STOP") != std::string::npos)
                    loc = NS;
                else { // Process constraint data
                    // Grab key and value
                    GetKeyValuePair(line, key, val);

                    // Set value for key
                    if (key.compare("T0") == 0)
                        t0 = std::stod(val);
                    else if (key.compare("A0") == 0) {
                        if (val.find("UNCONSTRAINED") != std::string::npos) {
                            a0 = 0.0;
                            a0_con_on = false;
                        }
                        else {
                            a0 = std::stod(val); 
                            a0_con_on = true;
                        }
                    }
                    else if (key.compare("E0") == 0) {
                        if (val.find("UNCONSTRAINED") != std::string::npos) {
                            e0 = 0.0;
                            e0_con_on = false;
                        }
                        else {
                            e0 = std::stod(val);
                            e0_con_on = true;
                        }
                    }
                    else if (key.compare("AOP0") == 0) {
                        if (val.find("UNCONSTRAINED") != std::string::npos) {
                            aop0 = 0.0;
                            aop0_con_on = false;
                        }
                        else {
                            aop0 = std::stod(val);
                            aop0_con_on = true;
                        }
                    }
                    else if (key.compare("TA0") == 0) {
                        if (val.find("UNCONSTRAINED") != std::string::npos) {
                            ta0 = 0.0;
                            ta0_con_on = false;
                        }
                        else {
                            ta0 = std::stod(val);
                            ta0_con_on = true;
                        }
                    }
                    else if (key.compare("ALPHA0") == 0) {
                        if (val.find("UNCONSTRAINED") != std::string::npos) {
                            alpha0 = 0.0;
                            alpha0_con_on = false;
                        }
                        else {
                            alpha0 = std::stod(val);
                            alpha0_con_on = true;
                        }
                    }
                    else if (key.compare("D_ALPHA0") == 0) {
                        if (val.find("UNCONSTRAINED") != std::string::npos) {
                            d_alpha0 = 0.0;
                            d_alpha0_con_on = false;
                        }
                        else {
                            d_alpha0 = std::stod(val);
                            d_alpha0_con_on = true;
                        }
                    }
                    else if (key.compare("L0") == 0) {
                        if (val.find("UNCONSTRAINED") != std::string::npos) {
                            L0 = 0.0;
                            L0_con_on = false;
                        }
                        else {
                            L0 = std::stod(val);
                            L0_con_on = true;
                        }
                    }
                    else if (key.compare("D_L0") == 0) {
                        if (val.find("UNCONSTRAINED") != std::string::npos) {
                            d_L0 = 0.0;
                            d_L0_con_on = false;
                        }
                        else {
                            d_L0 = std::stod(val);
                            d_L0_con_on = true;
                        }
                    }
                    else if (key.compare("RPF") == 0) {
                        if (val.find("UNCONSTRAINED") != std::string::npos) {
                            rpf = 0.0;
                            rpf_con_on = false;
                        }
                        else {
                            rpf = std::stod(val);
                            rpf_con_on = true;
                        }
                    }
                    else if (key.compare("TF_MIN") == 0) 
                        tfmin = std::stod(val);
                    else if (key.compare("TF_MAX") == 0)
                        tfmax = std::stod(val);
                }
            }
        }
    }
}