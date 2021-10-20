### A Pluto.jl notebook ###
# v0.16.3

using Markdown
using InteractiveUtils

# ╔═╡ f85cfc5e-26cc-11ec-3e25-7fc30ba989a3
begin
	using Pkg
	Pkg.activate()
	# uncomment the following cell to download this package
	# Pkg.add(url="https://github.com/MasonProtter/GeometricMatrixAlgebras.jl")
end

# ╔═╡ c45438a7-2d79-4293-8925-d6df8043f98b
using GeometricMatrixAlgebras

# ╔═╡ f656d603-394a-4d61-b401-7343b2831a2d
using GeometricMatrixAlgebras.Basis3D

# ╔═╡ b1b67006-4509-49d3-8cf1-655337e1fef2
using GeometricMatrixAlgebras.StaticArrays

# ╔═╡ 4b38ee8a-e21a-417d-9926-4814027dc3df
try
	using PlutoUI; TableOfContents()
catch e;
	nothing
end

# ╔═╡ 01483e51-f129-47bc-9b77-bc6bb9992790
md"# GeometricMatrixAlgebras.jl"

# ╔═╡ ea8e3d31-53e5-4315-8ec7-a71420102547
md"This package is a playground for learning about Geometric Algebra (GA), and intendes to leverage Julia's powerful Linear Algebra ecosystem to compute quntities in GA using matrices as a backend."

# ╔═╡ b95a10e2-9dbf-4f05-b234-940afe747f1b
md"GeometricMatrixAlgebras.jl exports two operators, `⋅` and `∧` (typed `\cdot<TAB>` and `\wedge<TAB>` respectively) for the inner and outer products."

# ╔═╡ dbf57e2d-f5ca-4d35-9715-f68212b9dd62
names(GeometricMatrixAlgebras) #the names exported by the package

# ╔═╡ c3113a5a-4c04-4410-bdc5-23063f1e8f5c
md"The `Basis3D` submodule defines everything we need to work with a 3D real geometric algebra."

# ╔═╡ 14527f6d-0141-4f07-9d51-a2f44fefdd97
names(Basis3D) #the names exported by the Basis3D submodule

# ╔═╡ c6cad1fc-d90e-4f7b-b775-8eb2930e35d5
md"### Basic algebraic manipulations"

# ╔═╡ 05d069f5-abd5-4ad9-b7fb-c36044440357
md"Basis vectors square to the identity and anti-commute"

# ╔═╡ 6943e3c4-e11e-4ac1-8048-b2d8ec21eb78
σ1 * σ1

# ╔═╡ 8fa351d1-43e5-49c3-b405-0348272ce47d
σ2 * σ1

# ╔═╡ 70ac954e-24ad-4ff9-ac2d-23cad48eff36
σ2 * (1 + σ2) * σ3

# ╔═╡ 0fd9ad2f-0bcb-47c2-ba0b-9ffed2bed383
md"Define a length 2 vector in the x-y plane"

# ╔═╡ 7dd18b97-5b6d-4588-bc71-d1a777c0e2f6
v = σ1 + σ2

# ╔═╡ ee28e702-554e-4d2b-8e2a-571dc3c142ff
md"Define a rotor which rotates a vector an angle of π/2 in the σ1σ2 plane"

# ╔═╡ 77087909-a612-41ef-a17c-50ad75b1e9d3
R = exp(π/2 * (σ1*σ2)/2)

# ╔═╡ 93fc24d2-b181-43d7-bea9-f43831116bd7
md"rotate v using R"

# ╔═╡ 92a6fd90-f33f-4c64-82a6-044d13092235
u = R'v*R

# ╔═╡ 82f5de38-335d-439b-85f8-1b3c20c4dce2
md"v and u are orthogonal"

# ╔═╡ ffdef128-cfbd-49c4-84f3-89d1681359f6
v ⋅ u

# ╔═╡ 1b9643ab-b651-49ff-931f-aa13d0e3e670
v ∧ u

# ╔═╡ e27363af-7d9a-4677-aa4a-6e3d7da929b6
md"## Turning a matrix into a MultiVector"

# ╔═╡ fd432b0b-7172-4c9b-9e92-e17652168a49
md"Lets say we have a function with which returns a basis:" 

# ╔═╡ cabb45d2-144d-415d-bed3-5e3100d40251
basis3d()

# ╔═╡ d8709e8f-6838-4db2-9fe1-311ff6a9d678
md"Here is how we project some random array onto our basis"

# ╔═╡ 26798fac-1d72-4f1b-b157-88d71ac4a66c
[1  5   9  13
 2  6  10  14
 3  7  11  15
 4  8  12  16] |> MultiVector{basis3d}

# ╔═╡ 98102281-ccfd-451d-bb73-24d0a9b53247
md"### Custom Bases"

# ╔═╡ ad0666de-061c-4ff3-a044-128ec2b55eea
md"To define your own custom basis, just make a function that returns all the elements of your basis as a named tuple. For example, here is how you could do that for the 2D geoemtric algebra"

# ╔═╡ 17d5fa39-2449-4072-ae0e-cd76b866c43b
function basis2d()
	𝟙 = [1 0 
		 0 1]
	e1 = [0 1
		  1 0]
	e2 = [1  0
		  0 -1]
	
	(;𝟙, e1, e2, e12 = e1*e2)
end

# ╔═╡ e4c32203-3f3e-452c-80f1-a3fa38b98226
e1 = MultiVector{basis2d}(basis2d().e1)

# ╔═╡ 17c54ad8-6e26-4e23-bd90-cdad1af13e50
e2 = MultiVector{basis2d}(basis2d().e2)

# ╔═╡ 5d2e93aa-ae46-4d91-931c-327b121f47b4
(e1 * e2)^2 

# ╔═╡ 0119a85d-0716-47b7-acea-707a1483d784
exp(e1*e2 * π)

# ╔═╡ 21f60dae-290b-4f52-99f4-6f18bc296d8e
md"In order to make this faster, we could use the StaticArrays.jl package and an `@generated function` to make sure the compiler does more work for us. (This is how `basis3d` was defined)"

# ╔═╡ be807a62-54e1-49e0-9813-e042774bbcfa
@generated function faster_basis2d()
	𝟙 = SA[1 0 
		   0 1]
	e1 = SA[0 1
		    1 0]
	e2 = SA[1  0
		    0 -1]
	quote
		(;𝟙=$𝟙, e1=$e1, e2=$e2, e12=$(e1*e2))
	end
end

# ╔═╡ 2d029832-b48d-420b-8bb6-f0f03c0d461a
md"CAVEAT: do not do the above trick with mutable arrays, it's only valid with static arrays."

# ╔═╡ 5ffe0cb2-dd38-4912-925e-5636c936108e
faster_basis2d()

# ╔═╡ 0c649a95-bc54-4fa6-844c-695869a7c45c
md"##### Contributions welcome!"

# ╔═╡ 0d22229d-7f90-4a87-8032-e7eac22f72e3
html"""
<style>
#launch_binder {
    display: none;
}
body.disable_ui main {
        max-width : 95%;
    }
@media screen and (min-width: 1081px) {
    body.disable_ui main {
        margin-left : 10px;
        max-width : 72%;
        align-self: flex-start;
    }
}
</style>
"""

# ╔═╡ Cell order:
# ╠═f85cfc5e-26cc-11ec-3e25-7fc30ba989a3
# ╟─4b38ee8a-e21a-417d-9926-4814027dc3df
# ╟─01483e51-f129-47bc-9b77-bc6bb9992790
# ╟─ea8e3d31-53e5-4315-8ec7-a71420102547
# ╠═c45438a7-2d79-4293-8925-d6df8043f98b
# ╟─b95a10e2-9dbf-4f05-b234-940afe747f1b
# ╠═dbf57e2d-f5ca-4d35-9715-f68212b9dd62
# ╠═f656d603-394a-4d61-b401-7343b2831a2d
# ╟─c3113a5a-4c04-4410-bdc5-23063f1e8f5c
# ╠═14527f6d-0141-4f07-9d51-a2f44fefdd97
# ╟─c6cad1fc-d90e-4f7b-b775-8eb2930e35d5
# ╟─05d069f5-abd5-4ad9-b7fb-c36044440357
# ╠═6943e3c4-e11e-4ac1-8048-b2d8ec21eb78
# ╠═8fa351d1-43e5-49c3-b405-0348272ce47d
# ╠═70ac954e-24ad-4ff9-ac2d-23cad48eff36
# ╟─0fd9ad2f-0bcb-47c2-ba0b-9ffed2bed383
# ╠═7dd18b97-5b6d-4588-bc71-d1a777c0e2f6
# ╟─ee28e702-554e-4d2b-8e2a-571dc3c142ff
# ╠═77087909-a612-41ef-a17c-50ad75b1e9d3
# ╟─93fc24d2-b181-43d7-bea9-f43831116bd7
# ╠═92a6fd90-f33f-4c64-82a6-044d13092235
# ╟─82f5de38-335d-439b-85f8-1b3c20c4dce2
# ╠═ffdef128-cfbd-49c4-84f3-89d1681359f6
# ╠═1b9643ab-b651-49ff-931f-aa13d0e3e670
# ╟─e27363af-7d9a-4677-aa4a-6e3d7da929b6
# ╟─fd432b0b-7172-4c9b-9e92-e17652168a49
# ╠═cabb45d2-144d-415d-bed3-5e3100d40251
# ╟─d8709e8f-6838-4db2-9fe1-311ff6a9d678
# ╠═26798fac-1d72-4f1b-b157-88d71ac4a66c
# ╟─98102281-ccfd-451d-bb73-24d0a9b53247
# ╟─ad0666de-061c-4ff3-a044-128ec2b55eea
# ╠═17d5fa39-2449-4072-ae0e-cd76b866c43b
# ╠═e4c32203-3f3e-452c-80f1-a3fa38b98226
# ╠═17c54ad8-6e26-4e23-bd90-cdad1af13e50
# ╠═5d2e93aa-ae46-4d91-931c-327b121f47b4
# ╠═0119a85d-0716-47b7-acea-707a1483d784
# ╟─21f60dae-290b-4f52-99f4-6f18bc296d8e
# ╠═b1b67006-4509-49d3-8cf1-655337e1fef2
# ╠═be807a62-54e1-49e0-9813-e042774bbcfa
# ╟─2d029832-b48d-420b-8bb6-f0f03c0d461a
# ╠═5ffe0cb2-dd38-4912-925e-5636c936108e
# ╟─0c649a95-bc54-4fa6-844c-695869a7c45c
# ╟─0d22229d-7f90-4a87-8032-e7eac22f72e3
