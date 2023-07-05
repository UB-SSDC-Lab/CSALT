
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