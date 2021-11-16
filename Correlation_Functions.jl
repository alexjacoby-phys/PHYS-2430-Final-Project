include("StateBuilder.jl")
include("Permutation_Matrices.jl")
include("InfoTools.jl")
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

for i in 1:9
    print(dot(psi,Bilin(5,i,9) *psi),"   ")
end
dot(psi,Bilin(1,5,9) *psi)
Bilin(1,5,9)









L=  16 #N^2
samples = 1000
paths = walks(L,samples)
pathmats = WalkMats(paths)
avglength = 0
for i in 1:samples
    avglength = avglength+sum(Distances(paths[i],4))*(1/samples)*(1/L)
end


avglength
α = 0.1
β = 1.2
γ = 0.#500/samples
psi = zeros(Complex{Float64},3^L)
for i in 1:samples
    dist = sum(Distances(paths[i],4))/L
    w = ((dist-β)/α)^2
    cw = ((avglength-β)/α)^2
    A = exp(cw-w+i*γ*im)
    #print(A,"     ")
    psi = psi + A*BuildVBS(pathmats[i,:,:])
    print("x")
end
normalize!(psi)

for i in 1:9
    print(dot(psi,Bilin(5,i,9) *psi),"   ")
end
dot(psi,Bilin(1,5,9) *psi)
Bilin(1,5,9)
