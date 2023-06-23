
#include "DebrisDeorbitPathObject.hpp"
#include "MessageInterface.hpp"

DebrisDeorbitPathObject::DebrisDeorbitPathObject() : UserPathFunction() {}

DebrisDeorbitPathObject::DebrisDeorbitPathObject(const DebrisDeorbitPathObject& copy) :
    UserPathFunction(copy)
{
    mu  = copy.mu;
    inc = copy.inc;
    m1  = copy.m1;
    m2  = copy.m2;
    mt  = copy.mt;
    ms  = copy.ms;
    l0  = copy.l0;
    k   = copy.k;
    c   = copy.c;
}

DebrisDeorbitPathObject& DebrisDeorbitPathObject::operator=(const DebrisDeorbitPathObject& copy)
{
    if (&copy == this)
        return *this;
   
    UserPathFunction::operator=(copy);
    mu  = copy.mu;
    inc = copy.inc;
    m1  = copy.m1;
    m2  = copy.m2;
    mt  = copy.mt;
    ms  = copy.ms;
    l0  = copy.l0;
    k   = copy.k;
    c   = copy.c;
   
   return *this;
}

DebrisDeorbitPathObject::~DebrisDeorbitPathObject() {}

void DebrisDeorbitPathObject::SetGravitationalParameter(Real muarg)
{
    mu = muarg;
}

void DebrisDeorbitPathObject::SetMaximumThrust(Real tMaxarg)
{
    tMax = tMaxarg;
}

void DebrisDeorbitPathObject::SetMassParameters(Real m1arg, Real m2arg, Real mtarg)
{
    // Set mass parameters
    m1  = m1arg;
    m2  = m2arg;
    mt  = mtarg;

    // Compute other mass parameters
    m   = m1 + m2 + mt;
    ms  = (m1 + 0.5*mt)*(m2 + 0.5*mt) / m - mt / 6.0;
}

void DebrisDeorbitPathObject::SetTetherParameters(Real l0arg, Real karg, Real carg)
{
    // Set parameters
    l0  = l0arg;
    k   = karg; 
    c   = carg;
}

void DebrisDeorbitPathObject::EvaluateFunctions()
{
    // Get state and control
    Rvector x = GetStateVector();
    Rvector u = GetControlVector();

    // Grab individual componants of state
    Real a      = x(0);
    Real e      = x(1);
    Real aop    = x(2);
    Real ta     = x(3);
    Real alpha  = x(4);
    Real dalpha = x(5);
    Real l      = x(6);
    Real dl     = x(7);

    // Grab individual componants of control
    Real Fx     = tMax*u(0)*u(1);
    Real Fz     = tMax*u(0)*u(2);

    // Compute accelerations due to thrust
    Real fu     = Fx / m;
    Real fr     = Fz / m;

    // Compute requirements
    Real p      = a*(1.0 - e*e);
    Real h      = sqrt(mu * p);
    Real R      = p / (1.0 + e*cos(ta));
    Real Om     = h / (R*R);
    Real dR     = R*e*Om*sin(ta) / (1.0 + e*cos(ta));
    Real dOm    = -2.0*Om / R * dR;
    Real l1     = (m2 + 0.5*mt) / m * l;

    // Generalized forces
    Real Qa     = (Fx*cos(alpha) - Fz*sin(alpha))*l1;
    Real Ql     = (Fx*sin(alpha) + Fz*cos(alpha))*l1 / l - k*(l - l0) - c*dl;

    // Equations of motion
    Rvector dx(8);
    Real sqrtMuP    = sqrt(mu * p);
    Real sqrtPoMu   = sqrt(p / mu);
    Real sinta      = sin(ta);
    Real costa      = cos(ta);
    Real cosa       = cos(alpha);
    Real invR2      = 1.0 / (R*R);
    Real invR3      = 1.0 / (R*R*R);
    Real invP       = 1.0 / p;
    Real invE       = 1.0 / e;
    Real OmpDalpha  = Om + dalpha;
    dx(0)           = 2.0*a*a/sqrtMuP*(e*sinta*fr + (1.0 + e*costa)*fu);
    dx(1)           = sqrtPoMu*(sinta*fr + ((1.0 + R*invP)*costa + e*R*invP)*fu);
    dx(2)           = invE * sqrtPoMu * (-costa*fr + (1.0 + R*invP)*sinta*fu);
    dx(3)           = sqrtMuP*invR2 + costa*invE*sqrtPoMu*fr - sinta*invE*(1.0 + R*invP)*sqrtPoMu*fu;
    dx(4)           = dalpha;
    dx(5)           = -dOm - 1.5*mu*invR3*sin(2.0*alpha) + Qa/(ms*l*l);
    dx(6)           = dl;
    dx(7)           = Ql/ms - (mu*(1.0 - 3.0*cosa*cosa)*invR3 - OmpDalpha*OmpDalpha)*l;

    // Set dynamics
    SetFunctions(DYNAMICS, dx);

    // Path constraints
    Rvector pathCon(2);
    Rvector pathConLB(2);
    Rvector pathConUB(2);
    pathCon(0)      = u(0);
    pathConLB(0)    = 0.0; 
    pathConUB(0)    = 1.0;
    pathCon(1)      = sqrt(u(1)*u(1) + u(2)*u(2));
    pathConLB(1)    = 1.0;
    pathConUB(1)    = 1.0;

    // Set path constraints
    SetFunctions(ALGEBRAIC, pathCon);
    SetFunctionBounds(ALGEBRAIC, LOWER, pathConLB);
    SetFunctionBounds(ALGEBRAIC, UPPER, pathConUB);
}

void DebrisDeorbitPathObject::EvaluateJacobians()
{
    if (true) {
    // Compute dynamics partials
    Rmatrix stateDynPartials    = ComputeStateDynamicsPartials();
    Rmatrix controlDynPartials  = ComputeControlDynamicsPartials();
    Rmatrix timeDynPartials     = ComputeTimeDynamicsPartials();

    // Set dynamics partials
    SetJacobian(DYNAMICS, STATE, stateDynPartials);
    SetJacobian(DYNAMICS, CONTROL, controlDynPartials);
    SetJacobian(DYNAMICS, TIME, timeDynPartials);

    // Compute path function partials
    Rmatrix stateConPartials    = ComputeStateConstraintPartials();
    Rmatrix controlConPartials  = ComputeControlConstraintPartials();
    Rmatrix timeConPartials     = ComputeTimeConstraintPartials();

    // Set control partials
    SetJacobian(ALGEBRAIC, STATE, stateConPartials);
    SetJacobian(ALGEBRAIC, CONTROL, controlConPartials);
    SetJacobian(ALGEBRAIC, TIME, timeConPartials);
    }
}

Rmatrix DebrisDeorbitPathObject::ComputeStateDynamicsPartials()
{
    // Get state and control
    Rvector x = GetStateVector();
    Rvector u = GetControlVector();

    // Grab individual componants of state
    Real a          = x(0);
    Real e          = x(1);
    Real aop        = x(2);
    Real ta         = x(3);
    Real alpha      = x(4);
    Real alpha_d    = x(5);
    Real L          = x(6);
    Real L_d        = x(7);

    // Grab individual componants of control
    Real u0         = u(0);
    Real u1         = u(1);
    Real u2         = u(2);

    Real t2 = cos(alpha);
    Real t3 = sin(alpha);
    Real t4 = cos(ta);
    Real t5 = sin(ta);
    Real t6 = a*a;
    Real t7 = alpha*2.0;
    Real t8 = e*e;
    Real t9 = e*e*e;
    Real t10 = ta*2.0;
    Real t21 = 1.0/a;
    Real t25 = 1.0/e;
    Real t27 = 1.0/m;
    Real t28 = 1.0/ms;
    Real t29 = 1.0/mu;
    Real t32 = mt/2.0;
    Real t11 = t2*t2;
    Real t12 = sin(t7);
    Real t13 = t4*t4;
    Real t14 = t4*t4*t4;
    Real t15 = sin(t10);
    Real t16 = t5*t5;
    Real t17 = t2*u1;
    Real t18 = e*t4;
    Real t19 = t3*u2;
    Real t20 = t4*u2;
    Real t22 = 1.0/t6;
    Real t23 = t21*t21*t21;
    Real t26 = 1.0/t8;
    Real t30 = e*t5*u2;
    Real t31 = t5*u1*2.0;
    Real t33 = t8-1.0;
    Real t38 = m2+t32;
    Real t51 = t5*t27*tMax*u0*u2;
    Real t24 = t22*t22;
    Real t34 = t11*3.0;
    Real t35 = -t19;
    Real t36 = -t31;
    Real t37 = t18+1.0;
    Real t41 = e*t13*u2*2.0;
    Real t42 = a*mu*t33;
    Real t45 = 1.0/t33;
    Real t49 = t8*t14*u2;
    Real t55 = a*t29*t33;
    Real t39 = t37*t37;
    Real t40 = t37*t37*t37;
    Real t43 = 1.0/t37;
    Real t46 = t45*t45;
    Real t47 = t45*t45*t45;
    Real t50 = t34-1.0;
    Real t52 = -t42;
    Real t53 = t17+t35;
    Real t56 = -t55;
    Real t44 = 1.0/t39;
    Real t48 = t46*t46;
    Real t54 = t43+1.0;
    Real t57 = sqrt(t52);
    Real t59 = sqrt(t56);
    Real t58 = 1.0/t57;
    Real t60 = 1.0/t59;
    Real t61 = t22*t39*t46*t57;
    Real t62 = alpha_d+t61;

    // Instantiate matrix for Jacobian
    Rmatrix A0(8,8);

    // Set entries of Jacobian to zero
    for (Integer i = 0; i < 8; i++)
        for (Integer j = 0; j < 8; j++)
            A0(i,j) = 0.0;

    // Set nonzero entries of Jacobian
    A0(0,0) = a*t27*t58*tMax*u0*(t30+u1+t18*u1)*3.0;
    A0(0,1) = t6*t27*t45*t58*tMax*u0*(e*u1+t4*u1+t5*u2)*-2.0;
    A0(0,3) = e*t6*t27*t58*tMax*u0*(t20-t5*u1)*2.0;
    A0(1,0) = t29*t33*t60*(t51+t27*tMax*u0*u1*(e*t43+t4*t54))*(-1.0/2.0);
    A0(1,1) = -a*t27*t29*t44*t60*tMax*u0*(t30-u1+t8*u1*2.0+t13*u1+t18*u1*2.0+t5*t8*t20*2.0+t4*t9*u1+t8*t13*u1*2.0+t9*t14*u1+t5*t9*t13*u2);
    A0(1,3) = t27*t44*t59*tMax*u0*(t20+t36+t41+t49+t5*t8*u1-t5*t18*u1*2.0-t5*t8*t13*u1);
    A0(2,0) = (t25*t29*t33*t60*(t20*t27*tMax*u0-t5*t27*t54*tMax*u0*u1))/2.0;
    A0(2,1) = a*t26*t27*t29*t44*t60*tMax*u0*(t20+t36+t41+t49-t8*u1*(t5-t5*t5*t5)-e*t15*u1*2.0+(t9*t15*u1)/2.0);
    A0(2,3) = t25*t59*(t51+t4*t27*t54*tMax*u0*u1+e*t16*t27*t44*tMax*u0*u1);
    A0(3,0) = t23*t39*t46*t57*-2.0-(mu*t22*t39*t45*t58)/2.0-(t20*t25*t27*t29*t33*t60*tMax*u0)/2.0+(t5*t25*t27*t29*t33*t54*t60*tMax*u0*u1)/2.0;
    A0(3,1) = e*t22*t39*t47*t57*-4.0+t4*t22*t37*t46*t57*2.0-e*mu*t21*t39*t46*t58-t20*t26*t27*t59*tMax*u0-a*t20*t27*t29*t60*tMax*u0+t5*t26*t27*t54*t59*tMax*u0*u1+a*t5*t27*t29*t54*t60*tMax*u0*u1+t4*t5*t25*t27*t44*t59*tMax*u0*u1;
    A0(3,3) = -t25*t51*t59-e*t5*t22*t37*t46*t57*2.0-t16*t27*t44*t59*tMax*u0*u1-t4*t25*t27*t54*t59*tMax*u0*u1;
    A0(4,5) = 1.0;
    A0(5,0) = mu*t24*t40*t47*(t12*3.0-e*t5*4.0)*(-3.0/2.0);
    A0(5,1) = (mu*t23*t39*t48*(t5*4.0-e*t12*1.8E+1+e*t15*8.0+t5*t8*2.0E+1-t4*t12*9.0+t9*t15*4.0-t4*t8*t12*9.0))/2.0;
    A0(5,3) = e*mu*t23*t39*t47*(e*-1.2E+1+t4*4.0+e*t13*1.6E+1+t5*t12*9.0)*(-1.0/2.0);
    A0(5,4) = mu*t23*t40*t47*(t11*2.0-1.0)*3.0-(t27*t28*t38*tMax*u0*(t2*u2+t3*u1))/L;
    A0(5,6) = -1.0/(L*L)*t27*t28*t38*t53*tMax*u0;
    A0(6,7) = 1.0;
    A0(7,0) = L*(mu*t24*t40*t47*t50*3.0+mu*t22*t39*t45*t58*t62*3.0);
    A0(7,1) = L*(e*mu*t23*t40*t48*t50*6.0-mu*t4*t23*t39*t47*t50*3.0+mu*t21*t46*t58*t62*(e*3.0+t4*2.0+e*t13*2.0+t4*t8*4.0+t9*t13)*2.0);
    A0(7,3) = L*(e*mu*t5*t23*t39*t47*t50*3.0-e*t5*t22*t37*t46*t57*t62*4.0);
    A0(7,4) = t27*t28*t38*t53*tMax*u0+L*mu*t2*t3*t23*t40*t47*6.0;
    A0(7,5) = L*t62*2.0;
    A0(7,6) = -k*t28+t62*t62-mu*t23*t40*t47*t50;
    A0(7,7) = -c*t28;
    return A0;
}

Rmatrix DebrisDeorbitPathObject::ComputeControlDynamicsPartials()
{
    // Get state and control
    Rvector x = GetStateVector();
    Rvector u = GetControlVector();

    // Grab individual componants of state
    Real a          = x(0);
    Real e          = x(1);
    Real aop        = x(2);
    Real ta         = x(3);
    Real alpha      = x(4);
    Real alpha_d    = x(5);
    Real L          = x(6);
    Real L_d        = x(7);

    // Grab individual componants of control
    Real u0         = u(0);
    Real u1         = u(1);
    Real u2         = u(2);

    Real t2 = cos(alpha);
    Real t3 = sin(alpha);
    Real t4 = cos(ta);
    Real t5 = sin(ta);
    Real t6 = a*a;
    Real t7 = e*e;
    Real t9 = 1.0/L;
    Real t10 = 1.0/e;
    Real t11 = 1.0/m;
    Real t12 = 1.0/ms;
    Real t13 = 1.0/mu;
    Real t14 = mt/2.0;
    Real t8 = e*t4;
    Real t15 = t7-1.0;
    Real t17 = m2+t14;
    Real t16 = t8+1.0;
    Real t18 = a*mu*t15;
    Real t22 = a*t13*t15;
    Real t19 = 1.0/t16;
    Real t20 = -t18;
    Real t23 = -t22;
    Real t21 = t19+1.0;
    Real t24 = 1.0/sqrt(t20);
    Real t25 = sqrt(t23);
    Real t26 = t4*t10*t11*t25*tMax*u0;
    Real t27 = t5*t10*t11*t21*t25*tMax*u0;

    // Instantiate matrix for jacobian
    Rmatrix A0(8,3);

    // Set entries of jacobian to zero
    for (Integer i = 0; i < 8; i++)
        for (Integer j = 0; j < 3; j++)
            A0(i,j) = 0.0;

    // Set nonzero jacobian entries
    A0(0,0) = t6*t11*t24*tMax*(u1+t8*u1+e*t5*u2)*2.0;
    A0(0,1) = t6*t11*t16*t24*tMax*u0*2.0;
    A0(0,2) = e*t5*t6*t11*t24*tMax*u0*2.0;
    A0(1,0) = t25*(t5*t11*tMax*u2+t11*tMax*u1*(e*t19+t4*t21));
    A0(1,1) = t11*t19*t25*tMax*u0*(e+t4*2.0+t4*t8);
    A0(1,2) = t5*t11*t25*tMax*u0;
    A0(2,0) = -t10*t25*(t4*t11*tMax*u2-t5*t11*t21*tMax*u1);
    A0(2,1) = t27;
    A0(2,2) = -t26;
    A0(3,0) = t4*t10*t11*t25*tMax*u2-t5*t10*t11*t21*t25*tMax*u1;
    A0(3,1) = -t27;
    A0(3,2) = t26;
    A0(5,0) = t9*t11*t12*t17*tMax*(t2*u1-t3*u2);
    A0(5,1) = t2*t9*t11*t12*t17*tMax*u0;
    A0(5,2) = -t3*t9*t11*t12*t17*tMax*u0;
    A0(7,0) = t11*t12*t17*tMax*(t2*u2+t3*u1);
    A0(7,1) = t3*t11*t12*t17*tMax*u0;
    A0(7,2) = t2*t11*t12*t17*tMax*u0;
    return A0;
}

Rmatrix DebrisDeorbitPathObject::ComputeTimeDynamicsPartials()
{
    Rmatrix A0(8,1);
    for (Integer i = 0; i < 8; i++)
        A0(i,0) = 0.0;
    return A0;
}

Rmatrix DebrisDeorbitPathObject::ComputeStateConstraintPartials()
{
    Rmatrix A0(2,8);
    for (Integer i = 0; i < 2; i++)
        for (Integer j = 0; j < 8; j++)
            A0(i,j) = 0.0;
    return A0;
}

Rmatrix DebrisDeorbitPathObject::ComputeControlConstraintPartials()
{
    // Get control
    Rvector u = GetControlVector();

    // Grab individual componants of control
    Real u0         = u(0);
    Real u1         = u(1);
    Real u2         = u(2);

    // Instantiate Jacobian
    Rmatrix A0(2,3);

    // First constraint
    A0(0,0)     = 1.0;
    A0(0,1)     = 0.0;
    A0(0,2)     = 0.0;

    // Second constraint
    Real ndir   = sqrt(u1*u1 + u2*u2);
    Real indir  = 1.0 / ndir;
    A0(1,0)     = 0.0;
    A0(1,1)     = u1 * indir;
    A0(1,2)     = u2 * indir;
    return A0;
}

Rmatrix DebrisDeorbitPathObject::ComputeTimeConstraintPartials()
{
    Rmatrix A0(2,1);
    A0(0,0) = 0.0;
    A0(1,0) = 0.0;
    return A0;
}