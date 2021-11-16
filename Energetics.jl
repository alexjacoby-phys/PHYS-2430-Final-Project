include("StateBuilder.jl")
include("Permutation_Matrices.jl")
include("Ops.jl")
include("ED.jl")
using Plots




L=  9 #N^2
samples = 10000
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
thet = []
compary = []
k=100
for i in 1:k
    θ = (2*pi)*(1-(i/k))
    gse = eigs(H(θ),nev = 1,which = :SR, maxiter = 1000)[1][1]
    energy = dot(psi,H(θ)*psi)
    comp = abs(real(energy-gse))
    push!(thet, θ)
    push!(compary,comp)
    print(comp, "          ")
end

eigs(H(pi/4),nev = 1,which = :SR, maxiter = 1000)[1]
plot(thet,compary)

thet
min(compary...)
