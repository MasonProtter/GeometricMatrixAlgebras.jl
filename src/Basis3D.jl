module Basis3D

using ..GeometricMatrixAlgebras: MultiVector, SA

export basis3d, 𝟙, σ1, σ2, σ3, σ23, σ31, σ12, σ123

@generated function basis3d()
    
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

const 𝟙  = MultiVector{basis3d}(basis3d().𝟙)

const σ1 = MultiVector{basis3d}(basis3d().σ1)
const σ2 = MultiVector{basis3d}(basis3d().σ2)
const σ3 = MultiVector{basis3d}(basis3d().σ3)

const σ23 = MultiVector{basis3d}(basis3d().σ23)
const σ31 = MultiVector{basis3d}(basis3d().σ31)
const σ12 = MultiVector{basis3d}(basis3d().σ12)

const σ123 = MultiVector{basis3d}(basis3d().σ123)

end
