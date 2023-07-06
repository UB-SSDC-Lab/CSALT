
# Define equations of motion function for easy differentiation
function ad_averaged_reduced_state_dynamics(x, u, t, μ, tMax, Isp, g0, τs, ws)
    dxtuple = averaged_reduced_state_dynamics(x..., t, u..., μ, tMax, Isp, g0, τs, ws)
    return SVector(dxtuple...)
end

function ad_averaged_reduced_state_dynamics_state_partial(
    p, f, g, h, k, m, t, δ, α1, α2, α3, μ, tMax, Isp, g0, τs, ws)

    # Construct state vector
    x = SVector(p, f, g, h, k, m)

    # Construct control vector
    u = SVector(δ, α1, α2, α3)

    # Compute partial
    ddxdx = ForwardDiff.jacobian(x -> ad_averaged_reduced_state_dynamics(x,u,t,μ,tMax,Isp,g0,τs,ws), x)

    # Place in single tuple and return
    return Tuple(ddxdx)
end

function ad_averaged_reduced_state_dynamics_control_partial(
    p, f, g, h, k, m, t, δ, α1, α2, α3, μ, tMax, Isp, g0, τs, ws)

    # Construct state vector
    x = SVector(p, f, g, h, k, m)

    # Construct control vector
    u = SVector(δ, α1, α2, α3)

    # Compute partial
    ddxdu = ForwardDiff.jacobian(u -> ad_averaged_reduced_state_dynamics(x,u,t,μ,tMax,Isp,g0,τs,ws), u)

    # Place in single tuple and return
    return Tuple(ddxdu)
end