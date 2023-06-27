#ifndef DEBRIS_DEORBIT_INITIAL_STATE_CONSTRAINT_HPP
#define DEBRIS_DEORBIT_INITIAL_STATE_CONSTRAINT_HPP

#include "OptimalControlFunction.hpp"

class DebrisDeorbitInitialStateConstraint : public OptimalControlFunction
{
public:
    // Default constructor
    DebrisDeorbitInitialStateConstraint(std::string funcName);

    // Copy constructor
    DebrisDeorbitInitialStateConstraint(const DebrisDeorbitInitialStateConstraint& copy);

    // Assignment operator
    DebrisDeorbitInitialStateConstraint& operator=(const DebrisDeorbitInitialStateConstraint& copy);

    // Default destructor
    virtual ~DebrisDeorbitInitialStateConstraint();

    // Evaluate optimal control functions
    virtual Rvector EvaluateFunctions();

    // Set constrained initial state
    void SetInitialStateConstraint(Rvector initialState);

    // Set constrained initial state idxs
    void SetInitialStateConstraintIdxs(IntegerArray idxsVec);

    // Set analytic jacobian on
    virtual void SetAnalyticJacobianOn(bool analyticJacobianOn);

    // Set upper and lower bounds
    virtual void SetLowerBounds(Rvector functionLB);
    virtual void SetUpperBounds(Rvector functionUB);

protected:
    // Evaluate analytic jacobians
    virtual void EvaluateAnalyticJacobian(VariableType varType, Integer pointIdx,
        bool &HasAnalyticJacobian, Rmatrix& jacArray);

    // Scale function bounds
    virtual void ScaleFunctionBounds();

    // Constrained initial states
    Rvector conInitialState;

    // Constrained initial state indecies
    IntegerArray conInitialStateIdxs;
};

#endif