include("StateBuilder.jl")
include("Permutation_Matrices.jl")
include("InfoTools.jl")




L=  9 #N^2
samples = 10^4
paths = walks(L,samples)
pathmats = WalkMats(paths)
avglength = 0
for i in 1:samples
    avglength = avglength+sum(Distances(paths[i],3))*(1/samples)*(1/L)
end

avglength
α = 0.15
β = 1.45
psi = zeros(Complex{Float64},3^L)
for i in 1:samples
    dist = sum(Distances(paths[i],3))/L
    w = (0.5*(dist-β)/α)^2
    cw = exp(-w)
    print(i, "     ")
    psi = psi + cw*BuildVBS(pathmats[i,:,:])
end
normalize!(psi)

DM = Sparse_DM_Build(psi)
eigvals(Matrix(Partial_Trace(DM, [1,2],9)))
eigvals(Matrix(Partial_Trace(DM, [1,4],9)))

eigvals(Matrix(Partial_Trace(DM, [1,3],9)))
eigvals(Matrix(Partial_Trace(DM, [1,7],9)))


eigvals(Matrix(Partial_Trace(DM, [1,9],9)))
