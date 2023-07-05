
#include "AveragedOrbitalElementsPathObject.hpp"
#include <iostream>

AveragedOrbitalElementsPathObject::AveragedOrbitalElementsPathObject() :
    n       (0),
    taus    (safe_eval("Vector{Float64}(undef, 0)")),
    ws      (safe_eval("Vector{Float64}(undef, 0)"))
{

}

AveragedOrbitalElementsPathObject::~AveragedOrbitalElementsPathObject()
{

}

void AveragedOrbitalElementsPathObject::Initialize(FunctionInputData *pd, PathFunctionContainer *pfc)
{
    // Ensure we have set gauss quadrature
    if (n == 0) {
        std::string errmsg = "ERROR ";
        errmsg += "initializing user path function in ";
        errmsg += "AveragedOrbitalElementsPathObject::Initialize. ";
        errmsg += "Cannot initialize path function before setting gauss ";
        errmsg += "quadrature.";
        throw LowThrustException(errmsg.c_str());
    }

    // Call UserPathFunction::Initialize()
    UserPathFunction::Initialize(pd, pfc);
}

void AveragedOrbitalElementsPathObject::SetGravitationalParameter(Real muarg)
{
    mu = muarg;
}

void AveragedOrbitalElementsPathObject::SetThrustParameters(Real tMaxarg, Real Isparg, Real g0arg)
{
    tMax    = tMaxarg;
    Isp     = Isparg;
    g0      = g0arg;
}

void AveragedOrbitalElementsPathObject::SetGaussQuadrature(Integer narg)
{
    // Set n member variable
    n = narg;

    // Get functions for getting nodes and weights
    auto get_taus   = Main["get_quadrature_nodes"];
    auto get_ws     = Main["get_quadrature_weights"];

    // Set weights and nodes
    taus    = get_taus(n);
    ws      = get_ws(n);
}

void AveragedOrbitalElementsPathObject::EvaluateFunctions()
{
    // Get states and controls
    Rvector states      = GetStateVector();
    Rvector controls    = GetControlVector();
    Real    t           = 0.0;

    // Grab individual states
    Real p  = states(0);
    Real f  = states(1);
    Real g  = states(2);
    Real h  = states(3);
    Real k  = states(4);
    Real m  = states(5);

    // Grab individual controls 
    Real sf = controls(0);
    Real a1 = controls(1);
    Real a2 = controls(2);
    Real a3 = controls(3);

    // Get EOM function
    auto eoms   = Main["averaged_reduced_state_dynamics"];

    // Call function
    //std::tuple<double, double, double, double, double, double> dx = eoms(p,f,g,h,k,m,t,sf,a1,a2,a3,mu,tMax,Isp,g0,taus,ws);
    auto dx = eoms(p,f,g,h,k,m,t,sf,a1,a2,a3,mu,tMax,Isp,g0,taus,ws);

    // Set dynamics 
    // Rvector dynamics(6, 
    //             std::get<0>(dx),
    //             std::get<1>(dx), 
    //             std::get<2>(dx), 
    //             std::get<3>(dx), 
    //             std::get<4>(dx), 
    //             std::get<5>(dx));
    Rvector dynamics(6);
    for (Integer i = 0; i < 6; i++)
        dynamics(i) = dx[i];

    SetFunctions(DYNAMICS, dynamics);

    // Path constraints
    Rvector pathCon(3);
    Rvector pathConLB(3);
    Rvector pathConUB(3);

    // Scale factor constraint
    pathCon(0)      = sf;
    pathConLB(0)    = 0.0;
    pathConUB(0)    = 1.0;

    // Direction unit norm constraint
    pathCon(1)      = sqrt(a1*a1 + a2*a2 + a3*a3);
    pathConLB(1)    = 1.0;
    pathConUB(1)    = 1.0;

    // Eccentricity constraint (to avoid singularity in EOMs)
    pathCon(2)      = sqrt(f*f + g*g);
    pathConLB(2)    = 0.0;
    pathConUB(2)    = 1.0;

    // Set path constraints
    SetFunctions(ALGEBRAIC, pathCon);
    SetFunctionBounds(ALGEBRAIC, LOWER, pathConLB);
    SetFunctionBounds(ALGEBRAIC, UPPER, pathConUB);
}

void AveragedOrbitalElementsPathObject::EvaluateJacobians()
{

}
