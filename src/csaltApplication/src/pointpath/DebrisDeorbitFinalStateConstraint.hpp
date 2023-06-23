
#ifndef DEBRIS_DEORBIT_FINAL_STATE_CONSTRAINT_HPP
#define DEBRIS_DEORBIT_FINAL_STATE_CONSTRAINT_HPP

#include "OptimalControlFunction.hpp"

class DebrisDeorbitFinalStateConstraint : public OptimalControlFunction
{
public:
    // Default constructor
    DebrisDeorbitFinalStateConstraint(std::string funcName);

    // Copy constructor
    DebrisDeorbitFinalStateConstraint(const DebrisDeorbitFinalStateConstraint& copy);

    // Assignment operator
    DebrisDeorbitFinalStateConstraint& operator=(const DebrisDeorbitFinalStateConstraint& copy);

    // Default destructor
    virtual ~DebrisDeorbitFinalStateConstraint();

    // Evaluate optimal control functions
    virtual Rvector EvaluateFunctions();

    // Set final state constraint
    void SetFinalStateConstraint(Real rpfarg);

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

    // Final constrained orbital radius
    Real rpf;
};

#endif