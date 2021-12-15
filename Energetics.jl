include("StateBuilder.jl")
include("Permutation_Matrices.jl")
include("Ops.jl")
include("ED.jl")
using Plots

gse = real(eigs(H(0.32175),nev = 1,which = :SR, maxiter = 1000)[1])[1]
maxen = real(eigs(H(0.32175),nev = 1,which = :LR, maxiter = 1000)[1])[1]-gse

L=  9 #N^2
samples = 2*10^4
paths = walks(L,samples)
pathmats = WalkMats(paths)
avglength = 0
for i in 1:samples
    avglength = avglength+sum(Distances(paths[i],3))*(1/samples)*(1/L)
end

avglength
α = 0.15
β = 1.
psi = zeros(Complex{Float64},3^L)
for i in 1:samples
    dist = sum(Distances(paths[i],3))/L
    w = (0.5*(dist-β)/α)^2
    cw = exp(-w)
    print(i, "     ")
    psi = psi + cw*BuildVBS(pathmats[i,:,:])
end
normalize!(psi)

Ham = H(0.32175)
(dot(psi,Ham*psi)-gse)/maxen
