# Cache quadrature rules
const squarerules = [QuadratureRule{2, RefCube}(order) for order = 1:5]
const cuberules = [QuadratureRule{3, RefCube}(order) for order = 1:5]
const trirules = [QuadratureRule{2, RefTetrahedron}(order) for order = 1:5]

function get_gaussrule(::Type{Dim{2}}, ::Type{RefTetrahedron}, order::Int)
    if order <= 5
        return trirules[order]
    else
        return QuadratureRule(Dim{2}, RefTetrahedron(), order)
    end
end

function get_gaussrule(::Type{Dim{2}}, ::Type{RefCube}, order::Int)
    if order <= 5
        return squarerules[order]
    else
        return QuadratureRule(Dim{2}, RefCube(), order)
    end
end

function get_gaussrule(::Type{Dim{3}}, ::Type{RefCube}, order::Int)
    if order <= 5
        return cuberules[order]
    else
        return QuadratureRule(Dim{3}, RefCube(), order)
    end
end
