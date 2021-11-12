include("StateBuilder.jl")
include("Permutation_Matrices.jl")
include("InfoTools.jl")
include("Ops.jl")
L= 9#N^2
samples = 100

psi = zeros(Complex{Float64},3^L)
paths = walks(L,samples)
pathmats = WalkMats(paths)


for i in 1:samples
    dist = sum(Distances(paths[i],3))
    psi = psi + exp(0.2*dist+5)*BuildVBS(pathmats[i,:,:])
    print("x")
end
normalize!(psi)

(norm(Collapse(3,9,9,psi)))^2




DM = Sparse_DM_Build(sparse(psi))
Partial_Trace(DM, [1,2,3,4,6,7,8,9],9)
for i in 1:9
    print(dot(psi,Bilin(1,i,9) *psi))
end
dot(psi,Bilin(1,7,9) *psi)
Bilin(1,2,16)

DM = Sparse_DM_Build(normalize(Collapse(1,5,9,psi)))
Partial_Trace(DM,[1,2,3,4,6,7,8,9],9)






#=
DM = Sparse_DM_Build(normalize!(Collapse(1,5,9,psi)))
for i in 0:7
    DM = Site_Trace(DM,9-i,9-i)
end

DM
=#
