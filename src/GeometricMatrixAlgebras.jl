module GeometricMatrixAlgebras

using LinearAlgebra, StaticArrays

export MultiVector, Basis3D, ∧, ⋅, basis

struct MultiVector{basis_maker, T <: AbstractMatrix}
    M::T
    MultiVector{basis_maker}(M::T) where {basis_maker, T <: AbstractMatrix} = new{basis_maker, T}(M)
end

basis(::MultiVector{basis_maker}) where {basis_maker} = basis_maker()


include("Basis3D.jl")

function decomp(basis::NamedTuple{basisnames}, M::AbstractMatrix) where {basisnames}
    N = size(M, 1)
    @assert N == size(M, 2)

    map(Tuple(basis)) do xi
        elem = (1//N) * (M ⋅ xi)
    end |> NamedTuple{basisnames}
end 

function Base.show(io::IO, x::MultiVector)
    GAprint(io, basis(x), x.M)
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
        print(io, "𝟘")
    else
        print(io, s)
    end
end



LinearAlgebra.:(⋅)(a::MultiVector, b::MultiVector) = 1//2 * (a*b + b*a) #typed \cdot<TAB>.
               (∧)(a::MultiVector, b::MultiVector) = 1//2 * (a*b - b*a) #dont confuse this with ^, the power operator

Base.isequal(a::MultiVector{B}, b::MultiVector{B}) where {B} = isequal(a.M, b.M)
Base.:(==)(a::MultiVector{B}, b::MultiVector{B}) where {B} = a.M == b.M
Base.isapprox(a::MultiVector{B}, b::MultiVector{B}; kwargs...) where {B} = isapprox(a.M, b.M; kwargs...)

for op ∈ (:*, :/, :^)
    @eval begin
        Base.$op(a::MultiVector{B}, b::MultiVector{B}) where {B} = MultiVector{B}($op(a.M, b.M))
        Base.$op(a::MultiVector{B}, b::Number) where {B} = MultiVector{B}($op(a.M, b))
        Base.$op(a::Number, b::MultiVector{B}) where {B} = MultiVector{B}($op(a, b.M))
    end
end
for op ∈ (:+, :-)
    @eval begin
        Base.$op(a::MultiVector{B}, b::MultiVector{B}) where {B} = MultiVector{B}($op(a.M, b.M))
        Base.$op(a::MultiVector{B}, b::Number) where {B} = MultiVector{B}($op(a.M, b * one(a.M)))
        Base.$op(a::Number, b::MultiVector{B}) where {B} = MultiVector{B}($op(a * one(b.M), b.M))
    end
end
for f ∈ (:+, :-, :log, :exp, :sin, :cos, :tan, :adjoint, :conj, :zero, :one)
    @eval Base.$f(a::MultiVector{B}) where {B} = MultiVector{B}($f(a.M))
end
Base.conj!(a::AbstractArray{<:MultiVector}) = a .= conj.(a)

# Base.:(*)(a::MultiVector, b::MultiVector) = MultiVector{basis(a)}(a.M * b.M)
# Base.:(*)(a::Number, b::MultiVector) = MultiVector{basis(b)}(a * b.M)
# Base.:(*)(a::MultiVector, b::Number) = MultiVector{basis(a)}(a.M * b)

# Base.:(/)(a::MultiVector, b::MultiVector) = MultiVector{basis(a)}(a.M / b.M)
# Base.:(/)(a::Number, b::MultiVector) = MultiVector{basis(b)}(a / b.M)
# Base.:(/)(a::MultiVector, b::Number) = MultiVector{basis(a)}(a.M / b)

# Base.:(+)(a::MultiVector, b::MultiVector) = MultiVector(a.M + b.M, a.basis)
# Base.:(-)(a::MultiVector, b::MultiVector) = MultiVector(a.M - b.M, a.basis)
# Base.:(^)(a::MultiVector, b::MultiVector) = MultiVector(a.M ^ b.M, a.basis)

# Base.:(+)(a::MultiVector) = a
# Base.:(-)(a::MultiVector) = MultiVector(-a.M, a.basis)

# Base.:(^)(a::Number, b::MultiVector) = MultiVector(a ^ b.M, b.basis)
# Base.:(^)(a::MultiVector, b::Number) = MultiVector(a.M ^ b, a.basis)

# Base.log(a::MultiVector) = MultiVector(log(a.M), a.basis)
# Base.exp(a::MultiVector) = MultiVector(exp(a.M), a.basis)
# Base.sin(a::MultiVector) = MultiVector(sin(a.M), a.basis)
# Base.cos(a::MultiVector) = MultiVector(cos(a.M), a.basis)
# Base.tan(a::MultiVector) = MultiVector(tan(a.M), a.basis)
# Base.adjoint(a::MultiVector) = MultiVector(a.M', a.basis)



end # module
