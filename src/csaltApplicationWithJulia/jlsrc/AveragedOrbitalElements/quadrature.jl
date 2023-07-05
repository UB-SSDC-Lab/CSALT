
# Helper functions for setting gauss quadrature
function get_quadrature_nodes(n)
    τs, ws = gausslegendre(n)
    return τs
end

function get_quadrature_weights(n)
    τs, ws = gausslegendre(n)
    return ws
end