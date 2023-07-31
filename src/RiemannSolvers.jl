module RiemannSolvers

using StaticArrays

abstract type AbstractRiemannSolver end

include("hllc.jl")

end
