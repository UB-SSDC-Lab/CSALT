
#include "AveragedOrbitalElementsPathObject.hpp"
#include <iostream>

//#define DEBUG_JULIA_PATH_FUNCTION

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
    auto dx = eoms(p,f,g,h,k,m,t,sf,a1,a2,a3,mu,tMax,Isp,g0,taus,ws);

    // Set dynamics 
    Rvector dynamics(6);
    for (Integer i = 0; i < 6; i++)
        dynamics(i) = dx[i];

#ifdef DEBUG_JULIA_PATH_FUNCTION
    std::cout << "DEBUG: Julia dynamics." << std::endl;
    for (Integer i = 0; i < 6; i++)
        std::cout << dynamics(i) << " " << (double) dx[i] << std::endl; 
#endif

    SetFunctions(DYNAMICS, dynamics);

    // Path constraints
    Rvector pathCon(2);
    Rvector pathConLB(2);
    Rvector pathConUB(2);

    // Scale factor constraint
    pathCon(0)      = sf;
    pathConLB(0)    = 0.0;
    pathConUB(0)    = 1.0;

    // Direction unit norm constraint
    pathCon(1)      = GmatMathUtil::Sqrt(a1*a1 + a2*a2 + a3*a3);
    pathConLB(1)    = 1.0;
    pathConUB(1)    = 1.0;

    // Set path constraints
    SetFunctions(ALGEBRAIC, pathCon);
    SetFunctionBounds(ALGEBRAIC, LOWER, pathConLB);
    SetFunctionBounds(ALGEBRAIC, UPPER, pathConUB);
}

void AveragedOrbitalElementsPathObject::EvaluateJacobians()
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

    // ===== Dynamics Partials

    // Get functions
    auto state_jac      = Main["averaged_reduced_state_dynamics_state_partial"];
    auto control_jac    = Main["averaged_reduced_state_dynamics_control_partial"];

    // Call functions
    auto ddxdx  = state_jac(p,f,g,h,k,m,t,sf,a1,a2,a3,mu,tMax,Isp,g0,taus,ws);
    auto ddxdu  = control_jac(p,f,g,h,k,m,t,sf,a1,a2,a3,mu,tMax,Isp,g0,taus,ws);

    // Allocate matricies for dynamics Jacobians
    Rmatrix dynStateJac(6,6), dynControlJac(6,4), dynTimeJac(6,1);

    // Create integer var for storing 
    Integer idx;

    // Fill dynamics state Jacobian
    idx = 0;
    for (Integer col = 0; col < 6; col++)
        for (Integer row = 0; row < 6; row++)
            dynStateJac(row,col) = ddxdx[idx++];

    // Fill dynamics control Jacobian
    idx = 0;
    for (Integer col = 0; col < 4; col++)
        for (Integer row = 0; row < 6; row++)
            dynControlJac(row,col) = ddxdu[idx++];

    // Fill dynamics time Jacobian
    for (Integer row = 0; row < 6; row++)
        dynTimeJac(row,0) = 0.0;

    // Set Jacobians
    SetJacobian(DYNAMICS, STATE,    dynStateJac);
    SetJacobian(DYNAMICS, CONTROL,  dynControlJac);
    SetJacobian(DYNAMICS, TIME,     dynTimeJac);

    // ===== Algebraic Partials

    // Allocate matricies for algebraic Jacobians
    Rmatrix algStateJac(2,6), algControlJac(2,4), algTimeJac(2,1);

    // Fill state Jacobian
    for (Integer col = 0; col < 6; col++)
        for (Integer row = 0; row < 2; row++)
            algStateJac(row,col) = 0.0;

    // Fill time Jacobian
    for (Integer row = 0; row < 2; row++)
        algTimeJac(row,0) = 0.0;

    // Fill control Jacobian
    Real invLen         = 1.0 / GmatMathUtil::Sqrt(a1*a1 + a2*a2 + a3*a3);
    algControlJac(0,0)  = 1.0;
    algControlJac(0,1)  = 0.0;
    algControlJac(0,2)  = 0.0;
    algControlJac(0,3)  = 0.0;
    algControlJac(1,0)  = 0.0;
    algControlJac(1,1)  = a1 * invLen; 
    algControlJac(1,2)  = a2 * invLen;
    algControlJac(1,3)  = a3 * invLen;

    // Set Jacobians
    SetJacobian(ALGEBRAIC, STATE,   algStateJac);
    SetJacobian(ALGEBRAIC, CONTROL, algControlJac);
    SetJacobian(ALGEBRAIC, TIME,    algTimeJac);
}
