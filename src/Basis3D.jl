module Basis3D

using ..GeometricMatrixAlgebras: MultiVector, SA
export basis3d, 𝟙, σ1, σ2, σ3

const basis3d = let
    
    𝟙 = SA[1 0 0 0
           0 1 0 0
           0 0 1 0
           0 0 0 1];
    
    σ1 = SA[0 0 1 0
            0 0 0 1
            1 0 0 0
            0 1 0 0];
    
    σ2 = SA[0  0  0  1
            0  0 -1  0
            0 -1  0  0
            1  0  0  0];
    
    σ3 = SA[1 0  0  0
            0 1  0  0
            0 0 -1  0
            0 0  0 -1]
    
    (;𝟙, #the scalar element
     σ1, σ2, σ3, #the three independant vectors in 3D
     σ23 = σ2*σ3, σ31 = σ3*σ1, σ12 = σ1*σ2, #the three independant bivectors in 3D
     σ123 = σ1*σ2*σ3 # the pseudoscalar element for 3D
     )
end

const 𝟙  = MultiVector(basis3d.𝟙, basis3d)
const σ1 = MultiVector(basis3d.σ1, basis3d)
const σ2 = MultiVector(basis3d.σ2, basis3d)
const σ3 = MultiVector(basis3d.σ3, basis3d)

end
