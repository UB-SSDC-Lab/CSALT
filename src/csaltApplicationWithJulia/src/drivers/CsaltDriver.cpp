
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

void CsaltDriver::SetInputFile(const std::string& inFile)
{
    inputFile = inFile;
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

            // Process solution
            ProcessSolution(traj);

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

void CsaltDriver::GetKeyValuePair(std::string s, std::string& key, std::string& val)
{
    // Get key by splitting by delimiter
    std::string del = "=";
    Integer pos = s.find(del);
    if (pos != std::string::npos) {
        key = s.substr(0, pos);

        // Remove white space
        key.erase(
            std::remove_if(
                key.begin(), 
                key.end(), 
                ::isspace), 
            key.end());

        // Remove substring from original string
        s.erase(0, pos + del.length());
    }
    else { // Should throw error here but just printing message for now
        std::cout << "Key value string missing equality sign delimiter!\n";
    }

    // Remove white space from remaining string and set value
    val = s;
    val.erase(
        std::remove_if(
            val.begin(),
            val.end(),
            ::isspace),
        val.end());
}

Rvector CsaltDriver::StringToRvector(std::string s)
{
    // Strip white space and brackets
    s.erase(std::remove_if(s.begin(),s.end(),::isspace),s.end());
    s.erase(std::remove(s.begin(),s.end(),'['),s.end());
    s.erase(std::remove(s.begin(),s.end(),']'),s.end());

    // Allocate vector to push values to
    Integer pos = 0;
    std::string sval;
    std::string del = ",";
    std::vector<Real> vals;
    while ((pos = s.find(del)) != std::string::npos)
    {
        // Get value and erase from string
        sval = s.substr(0, pos);
        s.erase(0, pos + del.length());

        // Push double to vector
        vals.push_back(std::stod(sval));
    }

    // Push final number to vector
    vals.push_back(std::stod(s));

    // Fill Rvector and return
    Rvector out(vals.size());
    for (Integer i = 0; i < vals.size(); i++)
        out(i) = vals.at(i);
    return out;
}

IntegerArray CsaltDriver::StringToIntegerArray(std::string s)
{
    // Strip white space and brackets
    s.erase(std::remove_if(s.begin(),s.end(),::isspace),s.end());
    s.erase(std::remove(s.begin(),s.end(),'['),s.end());
    s.erase(std::remove(s.begin(),s.end(),']'),s.end());

    // Allocate vector to push values to
    std::string sval;
    std::string del = ",";
    Integer pos = 0;
    IntegerArray out;
    while ((pos = s.find(del)) != std::string::npos)
    {
        // Get value  and erase from string
        sval = s.substr(0, pos);
        s.erase(0, pos + del.length());

        // Push double to vector
        out.push_back(std::stoi(sval));
    }

    // Push final number to vector
    out.push_back(std::stoi(s));

    return out;
}