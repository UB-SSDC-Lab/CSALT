
#include "AveragedOrbitalElementsDriver.hpp"
#include "AveragedOrbitalElementsPathObject.hpp"
#include "AveragedOrbitalElementsPointObject.hpp"
#include <string>

//#define DEBUG_JULIA_DRIVER

using namespace jluna;

AveragedOrbitalElementsDriver::AveragedOrbitalElementsDriver() : 
    CsaltDriver("AveragedOrbitalElements")
{
    // Set problem parameters
    SetParameters();
}

AveragedOrbitalElementsDriver::AveragedOrbitalElementsDriver(const std::string &testName, const Integer snoptConsoleOutputLevel) :
    CsaltDriver(testName, snoptConsoleOutputLevel)
{
    SetParameters();
}

AveragedOrbitalElementsDriver::~AveragedOrbitalElementsDriver()
{

}

void AveragedOrbitalElementsDriver::ProcessSolution(Trajectory* traj)
{
    // Get phase list
    std::vector<Phase*> phaseList = traj->GetPhaseList();

    // Get time vector and state/control arrays
    Rvector timeVector  = phaseList[0]->GetTimeVector();
    Rmatrix stateSol    = phaseList[0]->GetDecisionVector()->GetStateArray();
    Rmatrix controlSol  = phaseList[0]->GetDecisionVector()->GetControlArray();

    // Instantiate variabels for storing number of rows and cols
    Integer r, c;

    // Create string for constructing Julia code
    std::string code;

    // Create and fill Julia time vector
    r = timeVector.GetSize();
    code = "return zeros(" + std::to_string(r) + ")";
    Vector<Float64> times = safe_eval(code.c_str());
    for (size_t i = 0; i < r; i++)
        times[i] = timeVector(i);

    // Create and fill Julia state vector 
    stateSol.GetSize(r, c);
    code = "return zeros(" + std::to_string(r*c) + ")";
    Vector<Float64> states_vec = safe_eval(code.c_str());

    size_t idx = 0;
    for (size_t row = 0; row < r; row++)
        for (size_t col = 0; col < c; col++)
            states_vec[idx++] = stateSol(row,col);

#ifdef DEBUG_JULIA_DRIVER
    std::cout << code << std::endl;
    idx = 0;
    for (size_t row = 0; row < r; row++)
        for (size_t col = 0; col < c; col++)
            std::cout << stateSol(row,col) - (double) states_vec[idx++] << std::endl;
#endif

    // Create and fill Julia control vector
    controlSol.GetSize(r, c);
    code = "return zeros(" + std::to_string(r*c) + ")";
    Vector<Float64> controls_vec = safe_eval(code.c_str());

    idx = 0;
    for (size_t row = 0; row < r; row++)
        for (size_t col = 0; col < c; col++)
            controls_vec[idx++] = controlSol(row,col);

    // Call Julia process solution function
    auto proc_sol = Main["process_averaged_orbital_elements_solution"];
    proc_sol(times, states_vec, controls_vec);
}

void AveragedOrbitalElementsDriver::SetParameters()
{
    // Set file paths
    controlHistoryFile = "./../data/csalt/aao_sol.och";

    // Define problem parameters
    LU      = 6378.0;
    mu      = 3.986e5;
    TU      = GmatMathUtil::Sqrt(LU*LU*LU / mu);
    m0      = 1200.0; 
    P       = 5.0e3;
    Isp     = 1800.0;
    g0      = 9.80664;
    eta     = 0.55;
    tMax    = 2.0*eta*P / (g0 * Isp);

    // Set number of points use in Gauss quadrature
    n       = 10;

    // Set initial state vector constraint
    x0_con.SetSize(6);
    x0_con[0] = 11359.07 / LU;
    x0_con[1] = 0.7306;
    x0_con[2] = 0.0;
    x0_con[3] = 0.2539676;
    x0_con[4] = 0.0;
    x0_con[5] = m0;

    // Set final state vector constraint
    xf_con.SetSize(5);
    xf_con[0] = 42165.0 / LU;
    xf_con[1] = 0.01;
    xf_con[2] = 0.0;
    xf_con[3] = 8.726646282124052e-5;
    xf_con[4] = 0.0;

    // Scale variables
    mu      *= TU*TU / (LU*LU*LU);
    Isp     *= 1.0 / TU; 
    g0      *= TU*TU / (LU*1000.0);
    tMax    *= TU*TU / (LU*1000.0);
}

void AveragedOrbitalElementsDriver::SetPointPathAndProperties()
{
    // Instantiate point and path objects
    pathObject  = new AveragedOrbitalElementsPathObject();
    pointObject = new AveragedOrbitalElementsPointObject();

    // Set path object properties
    dynamic_cast<AveragedOrbitalElementsPathObject*>(pathObject)->SetGaussQuadrature(n);
    dynamic_cast<AveragedOrbitalElementsPathObject*>(pathObject)->SetGravitationalParameter(mu);
    dynamic_cast<AveragedOrbitalElementsPathObject*>(pathObject)->SetThrustParameters(tMax, Isp, g0);

    // Set point object properties
    dynamic_cast<AveragedOrbitalElementsPointObject*>(pointObject)->SetInitialStateConstraint(x0_con);
    dynamic_cast<AveragedOrbitalElementsPointObject*>(pointObject)->SetFinalStateConstraint(xf_con);

    // Set properties
    maxMeshRefinementCount      = 20;
    majorIterationsLimits[0]    = 100;
    totalIterationsLimits[0]    = 300000;
    optimizationMode            = StringArray(1, "Minimize");

    majorOptimalityTolerances.SetSize(1);
    feasibilityTolerances.SetSize(1);
    majorOptimalityTolerances(0)    = 1e-4;
    feasibilityTolerances(0)        = 1e-5;
}

void AveragedOrbitalElementsDriver::SetupPhases()
{
    // Instantiate LGR phase
    RadauPhase *phase   = new RadauPhase();

    // Set guess mode 
    std::string initialGuessMode    = "LinearUnityControl";

    // Set bounds
    Rvector state_LB(6);
    Rvector state_UB(6);
    state_LB(0) = 10000.0 / LU;
    state_UB(0) = 60000.0 / LU;
    state_LB(1) = -0.1;
    state_UB(1) =  0.8;
    state_LB(2) = -0.3;
    state_UB(2) =  0.3;
    state_LB(3) = -10.0;
    state_UB(3) =  10.0;
    state_LB(4) = -10.0;
    state_UB(4) =  10.0;
    state_LB(5) = 0.9 * m0;
    state_UB(5) = m0;

    Rvector control_LB(4);
    Rvector control_UB(4);
    control_LB(0) = 0.0;
    control_UB(0) = 1.0;
    for (Integer i = 1; i < 4; i++) {
        control_LB(i) = -1.0;
        control_UB(i) = 1.0;
    }

    Real time_LB = 0.0;
    Real time_UB = INF;

    // Set mesh properties
    Integer N, M;
    N = 10;
    M = 3;
    Real step = 2.0 / (N - 1.0);
    Rvector meshIntervalFractions(N);
    IntegerArray meshIntervalNumPoints;
    for (Integer i = 0; i < N - 1; i++) {
        meshIntervalFractions[i] = -1.0 + step*i;
        meshIntervalNumPoints.push_back(M);
    }
    meshIntervalFractions[N - 1] = 1.0;

    // Set initial and final guess for state
    Rvector x0_guess(6), xf_guess(6);
    for (Integer i = 0; i < 6; i++) {
        x0_guess(i) = x0_con(i);
        if (i != 5)
            xf_guess(i) = xf_con(i);
    }
    xf_guess(5) = 0.9*x0_guess(5);

    // Set phase properties
    phase->SetRelativeErrorTol(1e-6);
    phase->SetMaxPolynomialDegreeIncrease(3);
    phase->SetInitialGuessMode(initialGuessMode);
    phase->SetNumStateVars(6);
    phase->SetNumControlVars(4);
    phase->SetMeshIntervalFractions(meshIntervalFractions);
    phase->SetMeshIntervalNumPoints(meshIntervalNumPoints);
    phase->SetStateLowerBound(state_LB);
    phase->SetStateUpperBound(state_UB);
    phase->SetStateInitialGuess(x0_guess);
    phase->SetStateFinalGuess(xf_guess);
    phase->SetControlLowerBound(control_LB);
    phase->SetControlUpperBound(control_UB);
    phase->SetTimeLowerBound(0.0);
    phase->SetTimeUpperBound(time_UB);
    phase->SetTimeInitialGuess(0.0);
    phase->SetTimeFinalGuess(10.0*86400.0/TU);
    phaseList.push_back(phase);
}
