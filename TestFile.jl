include("StateBuilder.jl")
include("Permutation_Matrices.jl")
include("InfoTools.jl")
include("Ops.jl")
L= 9#N^2
samples = 200

psi = zeros(Complex{Float64},3^L)
paths = walks(L,samples)
pathmats = WalkMats(paths)


avglength = 0
for i in 1:samples
    avglength = avglength+sum(Distances(paths[i],3))*(1/samples)
end
avglength/9
α = 0.5
β = 1.5
for i in 1:samples
    dist = sum(Distances(paths[i],3))/L
    psi = psi + exp((avglength-β)^2/α)*exp(-(dist-β)^2/α)*BuildVBS(pathmats[i,:,:])
    print("x")
end
normalize!(psi)

(norm(Collapse(1,9,9,psi)))^2




DM = Sparse_DM_Build(sparse(psi))
Partial_Trace(DM, [1,2,3,4,6,7,8,9],9)
for i in 1:9
    print(dot(psi,Bilin(5,i,9) *psi),"   ")
end
dot(psi,Bilin(5,9,9) *psi)
Bilin(1,5,16)

DM = Sparse_DM_Build(normalize(Collapse(1,5,9,psi)))
Partial_Trace(DM,[1,2,3,4,6,7,8,9],9)





#=
DM = Sparse_DM_Build(normalize!(Collapse(1,5,9,psi)))
for i in 0:7
    DM = Site_Trace(DM,9-i,9-i)
end

DM
=#
