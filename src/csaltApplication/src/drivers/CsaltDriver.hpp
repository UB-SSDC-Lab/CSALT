
#ifndef CSALT_DRIVER_HPP
#define CSALT_DRIVER_HPP

#include "csalt.hpp"

class CsaltDriver
{
public:
    CsaltDriver(const std::string &probName);
    CsaltDriver(const std::string &probName, const Integer snoptConsoleOutputLevel);
    CsaltDriver(const CsaltDriver& driver);
    CsaltDriver& operator=(const CsaltDriver& driver);
    virtual ~CsaltDriver();

    enum Verbosity
    {
        SILENT,
        BASIC,
        VERBOSE,
        VERBOSE_DEBUG
    };

    virtual Integer Run();
protected:
    // Required method that sets the path & point objects
    virtual void SetPointPathAndProperties() = 0;

    // Required method that sets the phase(s)
    virtual void SetupPhases() = 0;

    // Methods for getting max error
    virtual Real GetMaxError(const Rvector &vec);
    virtual Real GetMaxError(const Rmatrix &mat);

    // Name of problem
    std::string probName;

    // Verbosity level
    Verbosity verbosity;

    // The trajectory
    Trajectory *traj;

    // List of phases
    std::vector<Phase*> phaseList;
    
    // User path and point functions
    UserPathFunction *pathObject;
    UserPointFunction *pointObject;

    // Cost function bounds
    Real costLowerBound;
    Real costUpperBound;

    // Snopt settings
    Rvector majorOptimalityTolerances;
    Rvector feasibilityTolerances;
    IntegerArray majorIterationsLimits;
    IntegerArray totalIterationsLimits;
    Integer maxMeshRefinementCount;
    StringArray optimizationMode;
    Integer snoptConsoleOutputLevel;

    // Output options 
    bool generateOptimizationOutput;
    std::string optimizationOutputFile;
    bool writeControlHistory;
    std::string controlHistoryFile;

    // A useful value
    static const Real INF;
};

#endif