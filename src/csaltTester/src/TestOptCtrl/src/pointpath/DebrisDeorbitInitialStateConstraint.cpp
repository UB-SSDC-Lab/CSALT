#include "DebrisDeorbitInitialStateConstraint.hpp"

DebrisDeorbitInitialStateConstraint::DebrisDeorbitInitialStateConstraint(std::string funcName) :
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

DebrisDeorbitInitialStateConstraint::DebrisDeorbitInitialStateConstraint(const DebrisDeorbitInitialStateConstraint& copy) :
    OptimalControlFunction(copy),
    conInitialState(copy.conInitialState),
    conInitialStateIdxs(copy.conInitialStateIdxs)
{}

DebrisDeorbitInitialStateConstraint& DebrisDeorbitInitialStateConstraint::operator=(const DebrisDeorbitInitialStateConstraint& copy) 
{
    if (this != &copy)
    {
        OptimalControlFunction::operator=(copy);
        conInitialState = copy.conInitialState;
        conInitialStateIdxs = copy.conInitialStateIdxs;
    }
    return *this;
}

DebrisDeorbitInitialStateConstraint::~DebrisDeorbitInitialStateConstraint() {}

Rvector DebrisDeorbitInitialStateConstraint::EvaluateFunctions()
{
    // Initialize constraint vector 
    Rvector con(numFunctions);

    // Get full state vector
    Rvector x = stateData.at(0);

    // Compute constraint
    for (Integer i = 0; i < numFunctions; i++) {
        con(i) = x(conInitialStateIdxs.at(i)) - conInitialState(i);
    }

    return con;
}

void DebrisDeorbitInitialStateConstraint::SetInitialStateConstraint(Rvector initialState)
{
    conInitialState.SetSize(initialState.GetSize());
    for (Integer i = 0; i < initialState.GetSize(); i++)
        conInitialState(i) = initialState(i);
}

void DebrisDeorbitInitialStateConstraint::SetInitialStateConstraintIdxs(IntegerArray idxsVec)
{
    conInitialStateIdxs.clear();
    for (Integer i = 0; i < idxsVec.size(); i++)
        conInitialStateIdxs.push_back(idxsVec.at(i));
}

void DebrisDeorbitInitialStateConstraint::SetAnalyticJacobianOn(bool analyticJacobianOn)
{
    analyticStateJacMap.resize(numPoints);
    if (analyticJacobianOn)
        analyticStateJacMap.at(0) = true;
    else 
        analyticStateJacMap.at(0) = false;
}

void DebrisDeorbitInitialStateConstraint::SetLowerBounds(Rvector functionLB)
{
    lowerBounds.SetSize(numFunctions);
    for (Integer i = 0; i < functionLB.GetSize(); i++)
        lowerBounds(i) = functionLB(i);
}

void DebrisDeorbitInitialStateConstraint::SetUpperBounds(Rvector functionUB)
{
    upperBounds.SetSize(numFunctions);
    for (Integer i = 0; i < functionUB.GetSize(); i++)
        upperBounds(i) = functionUB(i);
}

void DebrisDeorbitInitialStateConstraint::EvaluateAnalyticJacobian(VariableType varType, Integer pointIdx,
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

            // Set Jacobian matrix entries
            for (Integer i = 0; i < numFunctions; i++)
                for (Integer j = 0; j < 8; j++) 
                    if (j == conInitialStateIdxs.at(i)) 
                        jacArray(i,j) = 1.0;
                    else
                        jacArray(i,j) = 0.0;
            
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

void DebrisDeorbitInitialStateConstraint::ScaleFunctionBounds()
{
    ValidateFunctionBounds();
}