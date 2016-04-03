# Cache quadrature rules
const squarerules = [GaussQuadrature(Dim{2}, RefCube(), order) for order = 1:5]
const cuberules = [GaussQuadrature(Dim{3}, RefCube(), order) for order = 1:5]
const trirules = [GaussQuadrature(Dim{2}, RefTetrahedron(), order) for order = 1:5]

function get_gaussrule(::Type{Dim{2}}, ::RefTetrahedron, order::Int)
    if order <= 5
        return trirules[order]
    else
        return GaussQuadrature(Dim{2}, RefTetrahedron(), order)
    end
end

function get_gaussrule(::Type{Dim{2}}, ::RefCube, order::Int)
    if order <= 5
        return squarerules[order]
    else
        return GaussQuadrature(Dim{2}, RefCube(), order)
    end
end

function get_gaussrule(::Type{Dim{3}}, ::RefCube, order::Int)
    if order <= 5
        return cuberules[order]
    else
        return GaussQuadrature(Dim{3}, RefCube(), order)
    end
end
