include("StateBuilder.jl")
include("Permutation_Matrices.jl")
L= 16#N^2
psi = zeros(Complex{Float64},3^L)
paths = walks(L,200)
pathmats = WalkMats(paths)

BuildVBS(pathmats[15,:,:])
for i in 1:200
    psi = psi + BuildVBS(pathmats[i,:,:])
    print("x")
end
normalize!(psi)

sparse(I,3,3)
dot(psi,psi)
