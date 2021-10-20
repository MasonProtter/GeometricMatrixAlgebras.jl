using GeometricMatrixAlgebras, Test
using GeometricMatrixAlgebras.Basis3D

@testset "rotations" begin
    R = exp(π/2 * (σ1*σ2)/2) # a rotor which rotates a vector an angle of π/2 in the σ1σ2 plane
    v = σ1 + σ2              # a vector pointing 45 degrees in the x-y plane

    @test R'v*R ≈ -σ1 + σ2
end

@testset "Vectors of basis elements" begin
    @test [7, 11]' * [σ1, σ2] == (7σ1 + 11σ2)
    @test [7 11] * [σ1 σ2]' == fill(7σ1 + 11σ2, 1, 1)
end 

# ...
