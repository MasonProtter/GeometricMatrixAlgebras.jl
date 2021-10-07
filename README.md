# GeometricMatrixAlgebras

This package is a playground for learning about Geometric Algebra (GA), and intends to leverage Julia's powerful 
Linear Algebra ecosystem to compute quantities in GA using matrices as a backend.

This package is currently unregistered. To add it to a julia environment, do
```julia
using Pkg
Pkg.add(url="https://github.com/MasonProtter/GeometricMatrixAlgebras.jl")
```

To see this readme as a Pluto notebook, click [here](https://masonprotter.github.io/GeometricMatrixAlgebras.jl/)
___________

```julia
julia> using GeometricMatrixAlgebras

julia> using GeometricMatrixAlgebras.Basis3D # exports the basis elements for a 3D geometric algebra: 𝟙, σ1, σ2, σ3, basis3d

julia> σ1 * σ1
1 𝟙


julia> σ2 * σ1 
-1 σ12


julia> σ3 * σ2 * σ1
-1 σ123


julia> σ1 * σ2 * σ3 * σ1
1 σ23


julia> σ2 * (𝟙 + σ2) * σ3
1 σ3 + 1 σ23

julia> v = σ1 + σ2 # a vector pointing 45 degrees in the x-y plane
1 σ1 + 1 σ2

julia> R = exp(π/2 * (σ1*σ2)/2) # a rotor which rotates a vector an angle of π/2 in the σ1σ2 plane
0.7071067811865476 𝟙 + 2.7755575615628914e-17 σ3 + 0.7071067811865475 σ12

julia> u = R'v*R # rotate v using R
-0.9999999999999998 σ1 + 1.0000000000000git rem002 σ2 + -5.551115123125783e-17 σ31
```

GeometricMatrixAlgebras.jl exports two operators, `⋅` and `∧` (typed `\cdot<TAB>` and `\wedge` respectively) for the inner and outer products.

```julia 
julia> v ⋅ u # u and v are orthogonal
3.885780586188048e-16 𝟙

julia> v ∧ u
2.0 σ12
```

You can also wrap any matrix in an appropriate `MultiVector` given a basis, e.g.
```julia
julia> basis3d
(𝟙 = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1], σ1 = [0 0 1 0; 0 0 0 1; 1 0 0 0; 0 1 0 0], σ2 = [0 0 0 1; 0 0 -1 0; 0 -1 0 0; 1 0 0 0], σ3 = [1 0 0 0; 0 1 0 0; 0 0 -1 0; 0 0 0 -1], σ23 = [0 0 0 -1; 0 0 1 0; 0 -1 0 0; 1 0 0 0], σ31 = [0 0 1 0; 0 0 0 1; -1 0 0 0; 0 -1 0 0], σ12 = [0 -1 0 0; 1 0 0 0; 0 0 0 1; 0 0 -1 0], σ123 = [0 -1 0 0; 1 0 0 0; 0 0 0 -1; 0 0 1 0])

julia> [1  5   9  13
        2  6  10  14
        3  7  11  15
        4  8  12  16] |> MultiVector(basis3d) 
17//2 𝟙 + 17//2 σ1 + -5 σ3 + -3//2 σ23 + 3 σ31 + -3//2 σ123
```

Contributions welcome!
