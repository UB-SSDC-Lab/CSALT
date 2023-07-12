
#include "AveragedOrbitalElementsPointObject.hpp"

AveragedOrbitalElementsPointObject::AveragedOrbitalElementsPointObject()
{
    t0_con = 0.0;
}

AveragedOrbitalElementsPointObject::AveragedOrbitalElementsPointObject(const AveragedOrbitalElementsPointObject& copy)
{
    t0_con = copy.t0_con;
    x0_con = copy.x0_con;
    xf_con = copy.xf_con;
}

AveragedOrbitalElementsPointObject& AveragedOrbitalElementsPointObject::operator=(const AveragedOrbitalElementsPointObject& copy)
{
    if (&copy == this)
        return *this;

    UserPointFunction::operator=(copy);
    t0_con = copy.t0_con;
    x0_con = copy.x0_con;
    xf_con = copy.xf_con;
    return *this;
}

AveragedOrbitalElementsPointObject::~AveragedOrbitalElementsPointObject() {}

void AveragedOrbitalElementsPointObject::SetInitialStateConstraint(Rvector x0in)
{
    Integer n = x0in.GetSize();
    if (n != 6)
        throw LowThrustException("Initial state constraint vector must be of length 6.");

    x0_con.SetSize(n);
    for (Integer i = 0; i < n; i++)
        x0_con[i] = x0in[i];
}

void AveragedOrbitalElementsPointObject::SetFinalStateConstraint(Rvector xfin)
{
    Integer n = xfin.GetSize();
    if (n != 5)
        throw LowThrustException("Final state constraint vector must be of length 5.");

    xf_con.SetSize(n);
    for (Integer i = 0; i < n; i++)
        xf_con[i] = xfin[i];
}

void AveragedOrbitalElementsPointObject::EvaluateFunctions()
{
    // Get times
    Real t0 = GetInitialTime(0);
    Real tf = GetFinalTime(0);

    // Get states
    Rvector x0 = GetInitialStateVector(0);
    Rvector xf = GetFinalStateVector(0);

    // Instantiate vectors for constraint and bounds
    Rvector algF(12);
    Rvector algF_L(12);
    Rvector algF_U(12);

    //Rvector algF(11);
    //Rvector algF_L(11);
    //Rvector algF_U(11);

    //Rvector algF(10);
    //Rvector algF_L(10);
    //Rvector algF_U(10);

    // Set initial time constraint
    algF[0]     = t0 - t0_con;
    algF_L[0]   = 0.0;
    algF_U[0]   = 0.0;

    // Set initial state constraints
    for (Integer i = 0; i < 6; i++) {
        algF[i + 1]     = x0[i] - x0_con[i];
        algF_L[i + 1]   = 0.0;
        algF_U[i + 1]   = 0.0;
    }

    // Set final state constraint
    for (Integer i = 0; i < 5; i++) {
        algF[i + 7]     = xf[i] - xf_con[i];
        algF_L[i + 7]   = 0.0;
        algF_U[i + 7]   = 0.0;
    }

    //for (Integer i = 0; i < 4; i++) {
    //    if (i == 0)
    //        algF[i + 7]     = xf[i] - xf_con[i];
    //    else 
    //        algF[i + 7]     = xf[i + 1] - xf_con[i + 1];
    //    algF_L[i + 7]       = 0.0;
    //    algF_U[i + 7]       = 0.0;
    //}

    //Real a_f    = xf[0] / (1.0 - xf[1]*xf[1] - xf[2]*xf[2]);
    //Real e_f    = GmatMathUtil::Sqrt(xf[1]*xf[1] + xf[2]*xf[2]);
    //Real i_f    = GmatMathUtil::ATan2(2.0*GmatMathUtil::Sqrt(xf[3]*xf[3] + xf[4]*xf[4]), 1.0 - xf[3]*xf[3] - xf[4]*xf[4]);
    //algF[7]     = a_f -  42165.0 / 6378.0;
    //algF[8]     = e_f -  0.01;
    //algF[9]     = i_f -  0.01 * GmatMathConstants::PI / 180.0;
    //for (Integer i = 0; i < 3; i++) {
    //    algF_L[i + 7]  = 0.0;
    //    algF_U[i + 7]  = 0.0;
    //}

    // Set function and bounds
    SetFunctions(ALGEBRAIC, algF);
    SetFunctionBounds(ALGEBRAIC, LOWER, algF_L);
    SetFunctionBounds(ALGEBRAIC, UPPER, algF_U);

    // Set const function
    Rvector costFunc(1, tf);
    SetFunctions(COST, costFunc);
}

void AveragedOrbitalElementsPointObject::EvaluateJacobians()
{

}