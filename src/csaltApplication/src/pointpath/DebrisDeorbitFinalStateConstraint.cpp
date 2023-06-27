
#include "DebrisDeorbitFinalStateConstraint.hpp"

DebrisDeorbitFinalStateConstraint::DebrisDeorbitFinalStateConstraint(std::string funcName) :
    OptimalControlFunction(funcName) 
{
    // Set number of points
    numPoints = 1;

    // Resize variables for state, control, time, and parameter dependency maps
    stateDepMap.resize(numPoints);
    controlDepMap.resize(numPoints);
    timeDepMap.resize(numPoints);
    paramDepMap.resize(numPoints);

    // Set dependencies
    stateDepMap.at(0)   = true;
    controlDepMap.at(0) = false;
    timeDepMap.at(0)    = false;
    paramDepMap.at(0)   = false;
}

DebrisDeorbitFinalStateConstraint::DebrisDeorbitFinalStateConstraint(const DebrisDeorbitFinalStateConstraint& copy) :
    OptimalControlFunction(copy),
    rpf(copy.rpf) 
{}

DebrisDeorbitFinalStateConstraint& DebrisDeorbitFinalStateConstraint::operator=(const DebrisDeorbitFinalStateConstraint& copy) 
{
    if (this != &copy)
    {
        OptimalControlFunction::operator=(copy);
        rpf = copy.rpf;
    }
    return *this;
}

DebrisDeorbitFinalStateConstraint::~DebrisDeorbitFinalStateConstraint() {}

Rvector DebrisDeorbitFinalStateConstraint::EvaluateFunctions()
{
    // Initialize constraint vector 
    Rvector con(1);

    // Get full state vector
    Rvector x = stateData.at(0);

    // Compute constraint
    Real rp = x(0)*(1.0 - x(1));
    con(0) = rp - rpf;

    return con;
}

void DebrisDeorbitFinalStateConstraint::SetFinalStateConstraint(Real rpfarg)
{
    rpf = rpfarg;
}

void DebrisDeorbitFinalStateConstraint::SetAnalyticJacobianOn(bool analyticJacobianOn)
{
    analyticStateJacMap.resize(numPoints);
    if (analyticJacobianOn)
        analyticStateJacMap.at(0) = true;
    else 
        analyticStateJacMap.at(0) = false;
}

void DebrisDeorbitFinalStateConstraint::SetLowerBounds(Rvector functionLB)
{
    lowerBounds.SetSize(numFunctions);
    for (Integer i = 0; i < functionLB.GetSize(); i++)
        lowerBounds(i) = functionLB(i);
}

void DebrisDeorbitFinalStateConstraint::SetUpperBounds(Rvector functionUB)
{
    upperBounds.SetSize(numFunctions);
    for (Integer i = 0; i < functionUB.GetSize(); i++)
        upperBounds(i) = functionUB(i);
}

void DebrisDeorbitFinalStateConstraint::EvaluateAnalyticJacobian(VariableType varType, Integer pointIdx,
    bool &hasAnalyticJacobian, Rmatrix& jacArray)
{
    ValidatePointIdx(pointIdx);
    hasAnalyticJacobian = false;

    switch (varType)
    {
        case STATE:
        {
            // Set analytical Jacobain flag
            hasAnalyticJacobian = true;

            // Set size of Jacobian matrix
            jacArray.SetSize(numFunctions, 8);

            // Get state
            Rvector x = stateData.at(0);

            // Set Jacobian matrix entries
            jacArray(0,0) = 1.0 - x(1); 
            jacArray(0,1) = -x(0);
            for (Integer i = 2; i < 8; i++)
                jacArray(0,i) = 0.0;
            
            break;
        }
        case CONTROL: 
        {
            break;
        }
        case TIME: 
        {
            break;
        }
        case STATIC: 
        {
            break;
        }
    }
}

void DebrisDeorbitFinalStateConstraint::ScaleFunctionBounds()
{
    ValidateFunctionBounds();
}