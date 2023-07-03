
# Averaged equations of motion with the reduced averaged state vector (no true longitude) 
function averaged_reduced_state_dynamics(p, f, g, h, k, m, t, δ, α1, α2, α3, μ, tMax, Isp, g0, τs, ws)
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

    # Return averaged dynamics
    return pin, fin, gin, hin, kin, dmdt
end

# averaged_reduced_state_dynamics_ptr() = @cfunction(averaged_reduced_state_dynamics,
#             Tuple{Float64,Float64,Float64,Float64,Float64,Float64},
#            (Float64,Float64,Float64,Float64,Float64,Float64,
#             Float64,Float64,Float64,Float64,Float64,Float64,
#             Float64,Float64,Float64)