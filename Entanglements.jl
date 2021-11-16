include("StateBuilder.jl")
include("Permutation_Matrices.jl")
include("InfoTools.jl")




L=  9 #N^2
samples = 100
paths = walks(L,samples)
pathmats = WalkMats(paths)
avglength = 0
for i in 1:samples
    avglength = avglength+sum(Distances(paths[i],3))*(1/samples)*(1/L)
end


avglength
α = 0.1
β = 1.2
γ = 500/samples
psi = zeros(Complex{Float64},3^L)
for i in 1:samples
    dist = sum(Distances(paths[i],3))/L
    w = ((dist-β)/α)^2
    cw = ((avglength-β)/α)^2
    A = exp(cw-w+i*γ*im)
    #print(A,"     ")
    psi = psi + A*BuildVBS(pathmats[i,:,:])
    print("x")
end
normalize!(psi)









DM = Sparse_DM_Build(psi)
eigen(Matrix(Partial_Trace(DM, [2,3,4,5,6,7,8,9],9)))
eigen(Matrix(Partial_Trace(DM, [1,3,4,5,6,7,8,9],9)))
eigen(Matrix(Partial_Trace(DM, [1,2,4,5,6,7,8,9],9)))
eigen(Matrix(Partial_Trace(DM, [1,2,3,5,6,7,8,9],9)))
eigen(Matrix(Partial_Trace(DM, [1,2,3,4,6,7,8,9],9)))
eigen(Matrix(Partial_Trace(DM, [1,2,3,4,5,7,8,9],9)))
eigen(Matrix(Partial_Trace(DM, [1,2,3,4,5,6,8,9],9)))
eigen(Matrix(Partial_Trace(DM, [1,2,3,4,5,6,7,9],9)))
eigen(Matrix(Partial_Trace(DM, [1,2,3,4,5,6,7,8],9)))

eigen(Matrix(Partial_Trace(DM, [4,5,6,7,8,9],9)))
eigen(Matrix(Partial_Trace(DM, [1,2,3,4,5,6],9)))






DM = Sparse_DM_Build(normalize(Collapse(1,5,9,psi)))
Partial_Trace(DM,[1,2,3,4,6,7,8,9],9)
#=
