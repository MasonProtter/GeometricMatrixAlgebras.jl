module GeometricMatrixAlgebras

using LinearAlgebra, StaticArrays

export MultiVector, Basis3D, ⋅̂, ∧

struct MultiVector{T1 <: AbstractMatrix, T2 <: AbstractMatrix, basisnames, D}
    M::T1
    basis::NamedTuple{basisnames, NTuple{D, T2}}
end
MultiVector(basis::NamedTuple) = M -> MultiVector(M, basis)

include("Basis3D.jl")

function decomp(basis::NamedTuple{basisnames}, M::AbstractMatrix) where {basisnames}
    N = size(M, 1)
    @assert N == size(M, 2)

    map(Tuple(basis)) do xi
        elem = (1//N) * (M ⋅ xi)
    end |> NamedTuple{basisnames}
end 

function Base.show(io::IO, x::MultiVector{T1, T2,  D, basisnames}) where {T1, T2, D, basisnames}
    GAprint(io, x.basis, x.M)
end

GAprint(basis::NamedTuple) = M -> GAprint(stdout, basis, M)
GAprint(basis::NamedTuple, M::AbstractMatrix) = GAprint(stdout, basis, M)
function GAprint(io::IO, basis::NamedTuple{basisnames}, M::AbstractMatrix) where {basisnames}
    elems = map(decomp(basis, M) |> Tuple) do elem
        if elem isa Rational && denominator(elem) == 1
            numerator(elem)
        else
            elem
        end
    end
    
    s = join(("$(elems[i]) $(basisnames[i])" for i in eachindex(elems) if elems[i] != 0), " + ")
    if s == ""
        print("𝟘")
    else
        println(s)
    end
end



(⋅̂)(a, b) = 1//2 * (a*b + b*a) #typed \cdot<TAB>\hat<TAB>. I didn't use ⋅ because it's already taken by LinearAlgebra.jl
(∧)(a, b) = 1//2 * (a*b - b*a) #dont confuse this with ^, the power operator

Base.:(*)(a::MultiVector, b::MultiVector) = MultiVector(a.M * b.M, a.basis)
Base.:(*)(a::Number, b::MultiVector) = MultiVector(a * b.M, b.basis)
Base.:(*)(a::MultiVector, b::Number) = MultiVector(a.M * b, a.basis)

Base.:(/)(a::MultiVector, b::MultiVector) = MultiVector(a.M / b.M, a.basis)
Base.:(/)(a::Number, b::MultiVector) = MultiVector(a / b.M, b.basis)
Base.:(/)(a::MultiVector, b::Number) = MultiVector(a.M / b, a.basis)

Base.:(+)(a::MultiVector, b::MultiVector) = MultiVector(a.M + b.M, a.basis)
Base.:(-)(a::MultiVector, b::MultiVector) = MultiVector(a.M - b.M, a.basis)
Base.:(^)(a::MultiVector, b::MultiVector) = MultiVector(a.M ^ b.M, a.basis)

Base.:(^)(a::Number, b::MultiVector) = MultiVector(a ^ b.M, b.basis)
Base.:(^)(a::MultiVector, b::Number) = MultiVector(a.M ^ b, a.basis)

Base.log(a::MultiVector) = MultiVector(log(a.M), a.basis)
Base.exp(a::MultiVector) = MultiVector(exp(a.M), a.basis)
Base.sin(a::MultiVector) = MultiVector(sin(a.M), a.basis)
Base.cos(a::MultiVector) = MultiVector(cos(a.M), a.basis)
Base.tan(a::MultiVector) = MultiVector(tan(a.M), a.basis)
Base.adjoint(a::MultiVector) = MultiVector(a.M', a.basis)

end # module
