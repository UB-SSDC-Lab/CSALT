
#include "CsaltDriver.hpp"
#include "BaseException.hpp"
#include <iostream>

const Real CsaltDriver::INF = std::numeric_limits<Real>::infinity();

CsaltDriver::CsaltDriver(const std::string& probName) :
    probName                    (probName),
    verbosity                   (VERBOSE),
    traj                        (NULL),
    pathObject                  (NULL),
    pointObject                 (NULL),
    costLowerBound              (-INF),
    costUpperBound              (INF),
    maxMeshRefinementCount      (0),
    generateOptimizationOutput  (true),
    writeControlHistory         (true),
    snoptConsoleOutputLevel     (1)
{
    optimizationOutputFile      = probName + "Data.txt";
    controlHistoryFile          = probName + ".och";
    majorOptimalityTolerances   = Rvector(1, 1.0E-4);
    majorIterationsLimits       = IntegerArray(1, 3000);
    totalIterationsLimits       = IntegerArray(1, 300000);
    feasibilityTolerances       = Rvector(1, 1.0E-6);
    optimizationMode            = StringArray(1, "Minimize");
}

CsaltDriver::CsaltDriver(const std::string& probName, const Integer snoptConsoleOutputLevel) :
    probName                    (probName),
    verbosity                   (VERBOSE),
    traj                        (NULL),
    pathObject                  (NULL),
    pointObject                 (NULL),
    costLowerBound              (-INF),
    costUpperBound              (INF),
    maxMeshRefinementCount      (0),
    generateOptimizationOutput  (true),
    writeControlHistory         (true),
    snoptConsoleOutputLevel     (snoptConsoleOutputLevel)
{
    optimizationOutputFile      = probName + "Data.txt";
    controlHistoryFile          = probName + ".och";
    majorOptimalityTolerances   = Rvector(1, 1.0E-4);
    majorIterationsLimits       = IntegerArray(1, 3000);
    totalIterationsLimits       = IntegerArray(1, 300000);
    feasibilityTolerances       = Rvector(1, 1.0E-6);
    optimizationMode            = StringArray(1, "Minimize");
}

CsaltDriver::~CsaltDriver()
{
    if (traj)
        delete traj;
    if (pathObject)
        delete pathObject;
    if (pointObject)
        delete pointObject;
    for (Integer ii = 0; ii < phaseList.size(); ii++)
        if (phaseList.at(ii))
            delete phaseList.at(ii);
}

CsaltDriver::CsaltDriver(const CsaltDriver& driver)
{
}

CsaltDriver& CsaltDriver::operator=(const CsaltDriver& driver)
{
    return *this;
}

Integer CsaltDriver::Run()
{
    // Instantiate return value
    Integer retval = 0;

    // Set output format
    std::string outFormat = "%16.9f ";

    // Configure message interface and log file
    ConsoleMessageReceiver *consoleMsg = ConsoleMessageReceiver::Instance();
    MessageInterface::SetMessageReceiver(consoleMsg);
    std::string outPath = "./";
    MessageInterface::SetLogFile(outPath + "CsaltLog.txt");

    // Set global format settings
    GmatGlobal *global = GmatGlobal::Instance();
    global->SetActualFormat(false, false, 16, 1, false);

    // Output message about problem
    if (verbosity != SILENT)
    {
        MessageInterface::ShowMessage("%s\n",
                GmatTimeUtil::FormatCurrentTime().c_str());
        MessageInterface::ShowMessage("\n*** Running the %s CSALT problem ***\n",
            probName.c_str());
    }

    // Try to perform optimization
    try 
    {
        // Create the trajectory
        traj    = new Trajectory();

        // Call point and path function setup method
        SetPointPathAndProperties();

        // Set point and path functions for trajectory
        if (pathObject)
            traj->SetUserPathFunction(pathObject);
        if (pointObject)
            traj->SetUserPointFunction(pointObject);

        // Set SNOPT options
        traj->SetMajorIterationsLimit(majorIterationsLimits);
        traj->SetTotalIterationsLimit(totalIterationsLimits);
        traj->SetOptimalityTolerances(majorOptimalityTolerances);
        traj->SetFeasibilityTolerances(feasibilityTolerances);
        traj->SetOptimizationMode(optimizationMode);
        traj->SetCostLowerBound(costLowerBound);
        traj->SetCostUpperBound(costUpperBound);

        // Set CSALT options
        traj->SetMaxMeshRefinementCount(maxMeshRefinementCount);
        traj->SetFailedMeshOptimizationAllowance(true);
        traj->SetMeshRefinementGuessMode("LastSolutionMostRecentMesh");
        traj->SetSnoptConsoleOutputLevel(this -> snoptConsoleOutputLevel);

        // Setup the phases
        SetupPhases();

        if (phaseList.size() == 0) // Must have at least one phase, quit run
            retval = -1;
        else 
        {
            // Set phase list for trajectory
            traj->SetPhaseList(phaseList);

            // Initialize the trajectory
            traj->Initialize();

            // Prepare for call to optimize
            Rvector dv      = traj->GetDecisionVector();
            Rvector C       = traj->GetCostConstraintFunctions();
            Rvector z       = dv;
            Rvector F(C.GetSize());
            Rvector xmul(dv.GetSize());
            Rvector Fmul(C.GetSize());
            Integer exitFlag;

            // Call optimize
            if (generateOptimizationOutput)
                traj->Optimize(z, F, xmul, Fmul, exitFlag, optimizationOutputFile);
            else 
                traj->Optimize(z, F, xmul, Fmul, exitFlag);

            // Write OCH file
            if (writeControlHistory)
                traj->WriteToFile(controlHistoryFile);

            if (verbosity != SILENT)
                MessageInterface::ShowMessage("*** END %s PROBLEM ***\n", probName.c_str());
        }
    }
    catch(BaseException &ex)
    {
        MessageInterface::ShowMessage("Caught a CSALT Exception:\n\n%s\n\n",
            ex.GetFullMessage().c_str());
    }

    return retval;
}

Real CsaltDriver::GetMaxError(const Rvector &vec)
{
    Real max = -INF;

    for (Integer ii = 0; ii < vec.GetSize(); ii++)
        if (vec(ii) > max)
            max = vec(ii);

    return max;
}

Real CsaltDriver::GetMaxError(const Rmatrix &mat)
{
    Real max = -INF;
    Integer r, c;
    mat.GetSize(r,c);

    for (Integer ii = 0; ii < r; ii++)
        for (Integer jj = 0; jj < c; jj++)
            if (mat(ii,jj) > max)
                max = mat(ii,jj);

    return max;
}