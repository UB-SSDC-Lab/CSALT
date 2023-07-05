
#include "AveragedOrbitalElementsPathObject.hpp"
#include <iostream>

AveragedOrbitalElementsPathObject::AveragedOrbitalElementsPathObject(Integer n) : n(n)
{

}

AveragedOrbitalElementsPathObject::~AveragedOrbitalElementsPathObject()
{

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

void AveragedOrbitalElementsPathObject::EvaluateFunctions()
{
    // Get states and controls
    Rvector states = GetStateVector();
    Rvector controls = GetControlVector();
    Real    t = 0.0;

    // Set states in Julia vector and zero derivatives
    //JLvector x(6), dx(6);
    //for (Integer i = 0; i < 6; i++) {
    //    x[i] = states[i];
    //    dx[i] = 0.0;
    //}

    // Get quadrature weights and nodes
    //JLvector taus(n), ws(n);
    //jl_function_t* set_gq = jl_get_function(jl_main_module, "fill_gausslegendre!");
    //jl_value_t* args_gq[] = {(jl_value_t*) taus.GetJuliaVector(),
    //                         (jl_value_t*) ws.GetJuliaVector(),
    //                         jl_box_int64(n)};
    //handle_jl_call(set_gq, args_gq, 3);

    // Grab individual controls 
    //Real sf = controls(0);
    //Real a1 = controls(1);
    //Real a2 = controls(2);
    //Real a3 = controls(3);

    // Compute dynamics with Julia function: averaged_reduced_state_dynamics!(dx, x, t, δ, α1, α2, α3, μ, tMax, Isp, g0, τs, ws)
    // Get pointer to function
    //jl_function_t* eoms  = jl_get_function(jl_main_module, "averaged_reduced_state_dynamics!");

    // Construct array of arguments
    //jl_value_t* args[]   = {(jl_value_t*) dx.GetJuliaVector(), 
    //                        (jl_value_t*) x.GetJuliaVector(), 
    //                        jl_box_float64(t),
    //                        jl_box_float64(sf), 
    //                        jl_box_float64(a1), 
    //                        jl_box_float64(a2), 
    //                        jl_box_float64(a3), 
    //                        jl_box_float64(mu),
    //                        jl_box_float64(tMax), 
    //                        jl_box_float64(Isp), 
    //                        jl_box_float64(g0), 
    //                        (jl_value_t*) taus.GetJuliaVector(), 
    //                        (jl_value_t*) ws.GetJuliaVector()};

    // Call Julia function with 13 arguments and check for exceptions
    //handle_jl_call(eoms, args, 13); 

    // Set dynamics 
    //Rvector dynamics(6);
    //for (Integer i = 0; i < 6; i++) {
    //    dynamics[i] = dx[i];
    //}
    //handle_julia_exception();
    //SetFunctions(DYNAMICS, dynamics);

    // Path constraints
    //Rvector pathCon(2);
    //Rvector pathConLB(2);
    //Rvector pathConUB(2);
    //pathCon(0)      = sf;
    //pathConLB(0)    = 0.0;
    //pathConUB(0)    = 1.0;
    //pathCon(1)      = sqrt(a1*a1 + a2*a2 + a3*a3);
    //pathConLB(1)    = 1.0;
    //pathConUB(1)    = 1.0;

    // Set path constraints
    //SetFunctions(ALGEBRAIC, pathCon);
    //SetFunctionBounds(ALGEBRAIC, LOWER, pathConLB);
    //SetFunctionBounds(ALGEBRAIC, UPPER, pathConUB);
}

void AveragedOrbitalElementsPathObject::EvaluateJacobians()
{

}
