using GeometricMatrixAlgebras, Test
using GeometricMatrixAlgebras.Basis3D

@testset "rotations" begin
    R = exp(π/2 * (σ1*σ2)/2) # a rotor which rotates a vector an angle of π/2 in the σ1σ2 plane
    v = σ1 + σ2              # a vector pointing 45 degrees in the x-y plane

    R'v*R ≈ -σ1 + σ2
end
