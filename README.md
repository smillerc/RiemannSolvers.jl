# RiemannSolvers.jl

[![Build Status](https://github.com/smil/RiemannSolvers.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/smil/RiemannSolvers.jl/actions/workflows/CI.yml?query=branch%3Amain)


`RiemannSolvers.jl` provides a variety of different approximate Riemann solvers for use in computational fluid dynamics problems. See [Wikipedia](https://en.wikipedia.org/wiki/Riemann_solver) or the book by Toro []


# Example

```julia

ρ = 1.0 # density
v⃗1d = (1.0,) # 1d velocity vector
v⃗2d = (1.0, 2.0) # 2d velocity vector
v⃗3d = (1.0, 2.0, 3.0) # 3d velocity vector
p = 2.0 # pressure
n̂1d = (-1.0,) # 1d edge normal vector
n̂2d = (1.0, 0.0) # 2d edge normal vector
n̂3d = (1.0, 0.0, 1.0) # 3d face normal vector

F1d = flux(HLLC(), ρ, v⃗1d, p, n̂1d)
F2d = flux(HLLC(), ρ, v⃗2d, p, n̂2d)
F3d = flux(HLLC(), ρ, v⃗3d, p, n̂3d)
```