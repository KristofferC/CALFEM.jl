__precompile__()

module CALFEM

using JuAFEM
using InplaceOps
using Devectorize

import Base: show

import JuAFEM: Square, Triangle, FunctionSpace, n_dim, ref_shape, value!, derivative!,
       inv_spec!, det_spec


# Elements
export spring1e, spring1s
export plani4e, plani8e, soli8e, plante
export plani4s, plani8s, soli8s, plants
export plani4f, plani8f, soli8f, plantf
export flw2i4e, flw2i8e, flw2te, flw3i8e
export flw2i4s, flw2i8s, flw2ts, flw3i8s
export bar2e, bar2s, bar2g

# Materials
export hooke, solveq, assem
export extract, coordxtr, topologyxtr
export statcon, gen_quad_mesh

include("materials/hooke.jl")
include("utilities/utilities.jl")
include("elements/elements.jl")

end # module
