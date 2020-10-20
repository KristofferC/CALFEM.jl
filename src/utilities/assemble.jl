using LinearAlgebra

"""
Assembles the the element stiffness matrix `Ke`
to the global stiffness matrix `K`.
"""
function assem(edof::Vector, K::AbstractMatrix, Ke::Matrix,f::Union{VecOrMat, Nothing}=nothing,fe::Union{Vector, Nothing}=nothing)
    (nr, nc) = size(Ke)
    if nr != nc
        throw(DimensionMismatch("Stiffness matrix is not square (#rows=$nr #cols=$nc)"))
    elseif length(edof) != nr
        len_edof = length(edof)
        throw(DimensionMismatch("Mismatch between sizes in edof and Ke (edof($len_edof) Ke($nr,$nc)"))
    end
    K[edof, edof] += Ke
    return K
end
