
# Averaged equations of motion with the reduced averaged state vector (no true longitude) 
function averaged_reduced_state_dynamics(p, f, g, h, k, m, t, δ, α1, α2, α3, μ, tMax, Isp, g0, τs, ws)
    # Put thrust directions in SVector
    α = SVector(α1,α2,α3)

    # Compute orbital period for current osculating orbit
    a   = p / (1.0 - f*f - g*g) 
    T0  = 2.0*π*sqrt(a*a*a / μ)

    # Approximate the averaging integral with quadrature
    L0  = -π
    Lf  = π
    pin = 0.0
    fin = 0.0
    gin = 0.0
    hin = 0.0
    kin = 0.0
    @inbounds for i in eachindex(τs)
        # Compute true longitude for current point in quadrature
        Li      = 0.5*(Lf - L0)*τs[i] + 0.5*(Lf + L0)

        # Compute MEE dynamics
        xi      = SVector(p,f,g,h,k,Li,m)
        dmee    = mee_dynamics(xi, t, δ, α, μ, tMax)

        # Compute s
        w       = 1.0 + f*cos(Li) + g*sin(Li)
        dLdtp   = sqrt(μ*p)*(w / p)^2
        s       = 1.0 / (T0 * dLdtp)

        # Compute integrand
        int     = 0.5*(Lf - L0)*s*dmee

        # Add summation
        pin    += ws[i]*int[1]
        fin    += ws[i]*int[2]
        gin    += ws[i]*int[3]
        hin    += ws[i]*int[4]
        kin    += ws[i]*int[5]
    end

    # Compute mass dynamics
    dmdt = -tMax*δ / (Isp * g0)

    # Set averaged dynamics
    return (pin, fin, gin, hin, kin, dmdt)
end

# MEE element set dynamics including only central body point mass force model
function mee_dynamics(x, t, δ, α, μ, tMax)
    # Grab states
    p,f,g,h,k,L,m = x

    # Compute control
    Δ       = (tMax / m)*δ*α

    # Compute requirements
    w       = 1.0 + f*cos(L) + g*sin(L)
    ss      = 1.0 + h*h + k*k
    κ       = h*sin(L) - k*cos(L)
    sqrtpmu = sqrt(p / μ)
    wInv    = 1.0 / w

    # Compute dynamics
    dpdt    = 2.0*p*wInv*sqrtpmu*Δ[2]
    dfdt    = sqrtpmu*sin(L)*Δ[1] + sqrtpmu*wInv*((w+1.0)*cos(L) + f)*Δ[2] - sqrtpmu*g*wInv*κ*Δ[3]
    dgdt    = -sqrtpmu*cos(L)*Δ[1] + sqrtpmu*wInv*((w+1.0)*sin(L) + g)*Δ[2] + sqrtpmu*f*wInv*κ*Δ[3]
    dhdt    = 0.5*sqrtpmu*ss*cos(L)*wInv*Δ[3]
    dkdt    = 0.5*sqrtpmu*ss*sin(L)*wInv*Δ[3]
    dLdt    = sqrtpmu*wInv*κ*Δ[3] + sqrt(μ*p)*w*w/(p*p)

    # Return dynamics
    return SVector(dpdt, dfdt, dgdt, dhdt, dkdt, dLdt)
end

function mee_dynamics_partials(x, t, δ, α, μ, tMax)
    # Compute partials
    ddpdx = dpdotdx(x, t, δ, α, μ, tMax)
    ddfdx = dfdotdx(x, t, δ, α, μ, tMax)
    ddgdx = dgdotdx(x, t, δ, α, μ, tMax)
    ddhdx = dhdotdx(x, t, δ, α, μ, tMax)
    ddkdx = dkdotdx(x, t, δ, α, μ, tMax)
    ddLdx = dLdotdx(x, t, δ, α, μ, tMax)

    # Place in full Jacobian matrix
    return SMatrix{6,7}(ddpdx[1], ddfdx[1], ddgdx[1], ddhdx[1], ddkdx[1], ddLdx[1],
                        ddpdx[2], ddfdx[2], ddgdx[2], ddhdx[2], ddkdx[2], ddLdx[2],
                        ddpdx[3], ddfdx[3], ddgdx[3], ddhdx[3], ddkdx[3], ddLdx[3],
                        ddpdx[4], ddfdx[4], ddgdx[4], ddhdx[4], ddkdx[4], ddLdx[4],
                        ddpdx[5], ddfdx[5], ddgdx[5], ddhdx[5], ddkdx[5], ddLdx[5],
                        ddpdx[6], ddfdx[6], ddgdx[6], ddhdx[6], ddkdx[6], ddLdx[6],
                        ddpdx[7], ddfdx[7], ddgdx[7], ddhdx[7], ddkdx[7], ddLdx[7])
end

function dpdotdx(x, t, delta, α, mu, tMax)
    # Grab states
    p,f,g,h,k,L,m = x

    # Grab controls 
    ur,ut,un = α

    t2 = cos(L);
    t3 = sin(L);
    t6 = 1.0/m;
    t7 = 1.0/mu;
    t4 = f*t2;
    t5 = g*t3;
    t8 = p*t7;
    t9 = t4+t5+1.0;
    t10 = sqrt(t8);
    t11 = 1.0/t9;
    t12 = t11*t11;

    return SVector(
        delta*t6*t10*t11*tMax*ut*3.0,
        delta*p*t2*t6*t10*t12*tMax*ut*-2.0,
        delta*p*t3*t6*t10*t12*tMax*ut*-2.0,
        0.0,
        0.0,
        delta*p*t6*t10*t12*tMax*ut*(f*t3-g*t2)*2.0,
        delta*p*(t6*t6)*t10*t11*tMax*ut*-2.0
    )
end

function dfdotdx(x, t, delta, α, mu, tMax)
    # Grab states
    p,f,g,h,k,L,m = x

    # Grab controls 
    ur,ut,un = α

    t2 = cos(L);
    t3 = sin(L);
    t10 = 1.0/m;
    t12 = 1.0/mu;
    t4 = f*t2;
    t5 = g*t2;
    t6 = k*t2;
    t7 = f*t3;
    t8 = g*t3;
    t9 = h*t3;
    t11 = t10*t10;
    t15 = p*t12;
    t13 = -t7;
    t14 = -t9;
    t16 = t4+t8+1.0;
    t18 = sqrt(t15);
    t17 = t16+1.0;
    t19 = t5+t13;
    t20 = t6+t14;
    t21 = 1.0/t18;
    t23 = 1.0/t16;
    t22 = t2*t17;
    t24 = t23*t23;
    t25 = f+t22;

    return SVector(
        (delta*t3*t10*t12*t21*tMax*ur)/2.0+(delta*t10*t12*t21*t23*t25*tMax*ut)/2.0+(delta*g*t10*t12*t20*t21*t23*tMax*un)/2.0,
        delta*t10*t18*t23*tMax*ut*(t2*t2+1.0)-delta*t5*t10*t18*t20*t24*tMax*un-delta*t2*t10*t18*t24*t25*tMax*ut,
        delta*t10*t18*t20*t23*tMax*un-delta*t8*t10*t18*t20*t24*tMax*un+delta*t2*t3*t10*t18*t23*tMax*ut-delta*t3*t10*t18*t24*t25*tMax*ut,
        -delta*t8*t10*t18*t23*tMax*un,
        delta*t5*t10*t18*t23*tMax*un,
        delta*t2*t10*t18*tMax*ur-delta*t10*t18*t23*tMax*ut*(t3*t17-t2*t19)-delta*g*t10*t18*t23*tMax*un*(h*t2+k*t3)-delta*t10*t18*t19*t24*t25*tMax*ut-delta*g*t10*t18*t19*t20*t24*tMax*un,
        -delta*t3*t11*t18*tMax*ur-delta*t11*t18*t23*t25*tMax*ut-delta*g*t11*t18*t20*t23*tMax*un
    )
end

function dgdotdx(x, t, delta, α, mu, tMax)
    # Grab states
    p,f,g,h,k,L,m = x

    # Grab controls 
    ur,ut,un = α

    t2 = cos(L);
    t3 = sin(L);
    t10 = 1.0/m;
    t12 = 1.0/mu;
    t4 = f*t2;
    t5 = g*t2;
    t6 = k*t2;
    t7 = f*t3;
    t8 = g*t3;
    t9 = h*t3;
    t11 = t10*t10;
    t15 = p*t12;
    t13 = -t7;
    t14 = -t9;
    t16 = t4+t8+1.0;
    t18 = sqrt(t15);
    t17 = t16+1.0;
    t19 = t5+t13;
    t20 = t6+t14;
    t21 = 1.0/t18;
    t23 = 1.0/t16;
    t22 = t3*t17;
    t24 = t23*t23;
    t25 = g+t22;

    return SVector(
        delta*t2*t10*t12*t21*tMax*ur*(-1.0/2.0)+(delta*t10*t12*t21*t23*t25*tMax*ut)/2.0-(delta*f*t10*t12*t20*t21*t23*tMax*un)/2.0,
        -delta*t10*t18*t20*t23*tMax*un+delta*t4*t10*t18*t20*t24*tMax*un+delta*t2*t3*t10*t18*t23*tMax*ut-delta*t2*t10*t18*t24*t25*tMax*ut,
        delta*t10*t18*t23*tMax*ut*(t3*t3+1.0)+delta*t7*t10*t18*t20*t24*tMax*un-delta*t3*t10*t18*t24*t25*tMax*ut,
        delta*t7*t10*t18*t23*tMax*un,
        -delta*t4*t10*t18*t23*tMax*un,
        delta*t3*t10*t18*tMax*ur+delta*t10*t18*t23*tMax*ut*(t2*t17+t3*t19)+delta*f*t10*t18*t23*tMax*un*(h*t2+k*t3)-delta*t10*t18*t19*t24*t25*tMax*ut+delta*f*t10*t18*t19*t20*t24*tMax*un,
        delta*t2*t11*t18*tMax*ur-delta*t11*t18*t23*t25*tMax*ut+delta*f*t11*t18*t20*t23*tMax*un
    )
end

function dhdotdx(x, t, delta, α, mu, tMax)
    # Grab states
    p,f,g,h,k,L,m = x

    # Grab controls 
    ur,ut,un = α

    t2 = cos(L);
    t3 = sin(L);
    t4 = h*h;
    t5 = k*k;
    t8 = 1.0/m;
    t9 = 1.0/mu;
    t6 = f*t2;
    t7 = g*t3;
    t10 = p*t9;
    t11 = t4+t5+1.0;
    t12 = t6+t7+1.0;
    t13 = sqrt(t10);
    t14 = 1.0/t12;
    t15 = t14*t14;

    return SVector(
        (delta*t2*t8*t9*t11*t14*tMax*un)/(t13*4.0),
        delta*(t2*t2)*t8*t11*t13*t15*tMax*un*(-1.0/2.0),
        delta*t2*t3*t8*t11*t13*t15*tMax*un*(-1.0/2.0),
        delta*h*t2*t8*t13*t14*tMax*un,
        delta*k*t2*t8*t13*t14*tMax*un,
        delta*t3*t8*t11*t13*t14*tMax*un*(-1.0/2.0)+(delta*t2*t8*t11*t13*t15*tMax*un*(f*t3-g*t2))/2.0,
        delta*t2*(t8*t8)*t11*t13*t14*tMax*un*(-1.0/2.0)
    )
end

function dkdotdx(x, t, delta, α, mu, tMax)
    # Grab states
    p,f,g,h,k,L,m = x

    # Grab controls 
    ur,ut,un = α

    t2 = cos(L);
    t3 = sin(L);
    t4 = h*h;
    t5 = k*k;
    t8 = 1.0/m;
    t9 = 1.0/mu;
    t6 = f*t2;
    t7 = g*t3;
    t10 = p*t9;
    t11 = t4+t5+1.0;
    t12 = t6+t7+1.0;
    t13 = sqrt(t10);
    t14 = 1.0/t12;
    t15 = t14*t14;

    return SVector(
        (delta*t3*t8*t9*t11*t14*tMax*un)/(t13*4.0),
        delta*t2*t3*t8*t11*t13*t15*tMax*un*(-1.0/2.0),
        delta*(t3*t3)*t8*t11*t13*t15*tMax*un*(-1.0/2.0),
        delta*h*t3*t8*t13*t14*tMax*un,
        delta*k*t3*t8*t13*t14*tMax*un,
        (delta*t2*t8*t11*t13*t14*tMax*un)/2.0+(delta*t3*t8*t11*t13*t15*tMax*un*(f*t3-g*t2))/2.0,
        delta*t3*(t8*t8)*t11*t13*t14*tMax*un*(-1.0/2.0)
    )
end

function dLdotdx(x, t, delta, α, mu, tMax)
    # Grab states
    p,f,g,h,k,L,m = x

    # Grab controls 
    ur,ut,un = α

    t2 = cos(L);
    t3 = sin(L);
    t4 = mu*p;
    t11 = 1.0/m;
    t12 = 1.0/mu;
    t13 = 1.0/(p*p);
    t5 = f*t2;
    t6 = g*t2;
    t7 = k*t2;
    t8 = f*t3;
    t9 = g*t3;
    t10 = h*t3;
    t16 = p*t12;
    t17 = sqrt(t4);
    t14 = -t8;
    t15 = -t10;
    t18 = t5+t9+1.0;
    t19 = sqrt(t16);
    t20 = t6+t14;
    t21 = t7+t15;
    t22 = t18*t18;
    t23 = 1.0/t18;
    t24 = 1.0/t22;

    return SVector(
        1.0/(p*p*p)*t17*t22*-2.0+(mu*t13*t22)/(t17*2.0)-(delta*t11*t12*t21*t23*tMax*un)/(t19*2.0),
        t2*t13*t17*t18*2.0+delta*t2*t11*t19*t21*t24*tMax*un,
        t3*t13*t17*t18*2.0+delta*t3*t11*t19*t21*t24*tMax*un,
        delta*t3*t11*t19*t23*tMax*un,
        -delta*t2*t11*t19*t23*tMax*un,
        t13*t17*t18*t20*2.0+delta*t11*t19*t23*tMax*un*(h*t2+k*t3)+delta*t11*t19*t20*t21*t24*tMax*un,
        delta*(t11*t11)*t19*t21*t23*tMax*un
    )
end

function s_partials(p, f, g, h, k, L, m, mu)
    t2 = cos(L);
    t3 = sin(L);
    t4 = mu*p;
    t5 = f*f;
    t6 = g*g;
    t7 = p*p;
    t8 = p*p*p;
    t11 = 1.0/3.141592653589793;
    t12 = 1.0/mu;
    t9 = f*t2;
    t10 = g*t3;
    t13 = 1.0/sqrt(t4);
    t14 = t5+t6-1.0;
    t15 = t9+t10+1.0;
    t16 = 1.0/t14;
    t17 = t16*t16*t16;
    t18 = 1.0/(t15*t15*t15);
    t19 = t8*t12*t17;
    t20 = -t19;
    t21 = 1.0/sqrt(t20);
    return SVector(
        0.0,
        (t7*t11*t13*t16*t18*t21*(f*3.0+t2*2.0+f*t10*3.0+t2*t5-t2*t6*2.0))/2.0,
        (t7*t11*t13*t16*t18*t21*(g*3.0+t3*2.0+g*t9*3.0-t3*t5*2.0+t3*t6))/2.0,
        0.0,
        0.0,
        t7*t11*t13*t18*t21*(f*t3-g*t2),
        0.0
    )
end