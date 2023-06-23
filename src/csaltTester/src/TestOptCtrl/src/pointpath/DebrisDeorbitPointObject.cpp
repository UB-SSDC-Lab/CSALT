
#include "DebrisDeorbitPointObject.hpp"
#include "MessageInterface.hpp"

DebrisDeorbitPointObject::DebrisDeorbitPointObject() : UserPointFunction() {}

DebrisDeorbitPointObject::DebrisDeorbitPointObject(const DebrisDeorbitPointObject& copy) :
    UserPointFunction(copy),
    t0(copy.t0),
    tfmin(copy.tfmin),
    tfmax(copy.tfmax)
{
}

DebrisDeorbitPointObject& DebrisDeorbitPointObject::operator=(const DebrisDeorbitPointObject& copy)
{
    if (&copy == this)
        return *this;

    UserPointFunction::operator=(copy);
    t0    = copy.t0;
    tfmin = copy.tfmin;
    tfmax = copy.tfmax;
    return *this;
}

DebrisDeorbitPointObject::~DebrisDeorbitPointObject() {}

void DebrisDeorbitPointObject::SetTimeConstraints(Real t0arg, Real tfminarg, Real tfmaxarg)
{
    t0      = t0arg;
    tfmin   = tfminarg;
    tfmax   = tfmaxarg;
}

void DebrisDeorbitPointObject::EvaluateFunctions()
{
    // Get times
    Real t0     = GetInitialTime(0);
    Real tf     = GetFinalTime(0);

    // Instantiate vector for constraints and bounds
    Rvector algF(2);
    Rvector algF_L(2);
    Rvector algF_U(2);
    
    // Time constraints
    algF(0)     = t0;
    algF(1)     = tf;

    // Set algebraic constraint bounds
    algF_L(0)   = t0;
    algF_U(0)   = t0;
    algF_L(1)   = tfmin;
    algF_U(1)   = tfmax;

    // Set algebraic constraints
    SetFunctions(ALGEBRAIC, algF);
    SetFunctionBounds(ALGEBRAIC, LOWER, algF_L);
    SetFunctionBounds(ALGEBRAIC, UPPER, algF_U);

    // Set cost function
    Rvector costFunc(1, tf);
    SetFunctions(COST, costFunc);
}

void DebrisDeorbitPointObject::EvaluateJacobians()
{
    // Do nothing here (use finite diff for Jac)
}