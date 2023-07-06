
function process_averaged_orbital_elements_solution(times, states_vec, controls_vec)
    # Construct path to folder to save data
    dir = joinpath(@__DIR__, "..", "..", "..", "..", 
            "application", "data", "csalt")

    # Construct state and control matricies
    states      = zeros(Int(length(states_vec) / 6), 6)
    controls    = zeros(Int(length(controls_vec) / 4), 4)
    
    # States
    idx = 1
    for col in axes(states, 2)
        for row in axes(states, 1)
            states[row,col] = states_vec[idx]
            idx += 1
        end
    end

    # Controls
    idx = 1
    for col in axes(controls, 2)
        for row in axes(controls, 1)
            controls[row,col] = controls_vec[idx]
            idx += 1
        end
    end

    # Write state and control history to file
    write_to_file(dir, times, states, controls)

    # Plot states
    plot_states(dir, times, states)

    # Plot controls
    plot_controls(dir, times, controls)

    return nothing
end

function write_to_file(dir, times, states, controls)
    # Construct state and time matrix
    sdata = hcat(times, states)

    # Write to file
    writedlm(joinpath(dir, "aoe_states.txt"), sdata, ',')

    # Construct control and time matrix
    cdata = hcat(times[1:end - 1], controls)

    # Write to file
    writedlm(joinpath(dir, "aoe_controls.txt"), cdata, ',')

    return nothing
end

function plot_states(dir, times, states)
    # Create figure for states
    fig, axs = subplots(6; sharex = true)
    fig.suptitle("Averaged Orbital Element States")

    # Plot states
    for i in 1:6
        # Plot line
        axs[i].plot(times, states[:,i])

        # Get y-axis label
        ylabel = @match i begin 
            1 => "p" 
            2 => "f"
            3 => "g"
            4 => "h"
            5 => "k"
            6 => "m"
        end

        # Set y-axis label
        axs[i].set_ylabel(ylabel)

        # Set x-axis label if i == 6
        if i == 6
            axs[i].set_xlabel("time")
        end
    end

    # Save figure
    savefig(joinpath(dir, "aoe_states.pdf"); dpi = 300)
end

function plot_controls(dir, times, controls)
    # Create figure for states
    fig, axs = subplots(4; sharex = true)
    fig.suptitle("Averaged Orbital Element Controls")

    # Plot states
    for i in 1:4
        # Plot line
        axs[i].plot(times[1:end - 1], controls[:,i])

        # Get y-axis label
        ylabel = @match i begin 
            1 => "sf" 
            2 => "a1"
            3 => "a2"
            4 => "a3"
        end

        # Set y-axis label
        axs[i].set_ylabel(ylabel)

        # Set x-axis label if i == 6
        if i == 4
            axs[i].set_xlabel("time")
        end
    end

    # Save figure
    savefig(joinpath(dir, "aoe_controls.pdf"); dpi = 300)
end