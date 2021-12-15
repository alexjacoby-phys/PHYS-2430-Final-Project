include("StateBuilder.jl")
include("Permutation_Matrices.jl")
include("InfoTools.jl")
include("Ops.jl")
include("ED.jl")
using Plots
using LaTeXStrings
using DelimitedFiles


corfunc = []
abscorfunc = []
dvecnew = []


for i in 1:2
    L=  9 #N^2
    samples = 10^5
    paths = walks(L,samples)
    pathmats = WalkMats(paths)
    avglength = 0
    for i in 1:samples
        avglength = avglength+sum(Distances(paths[i],3))*(1/samples)*(1/L)
    end

    avglength
    α = 0.1
    β = 1.2
    psi = zeros(Complex{Float64},3^L)
    for i in 1:samples
        dist = sum(Distances(paths[i],3))/L
        w = (0.5*(dist-β)/α)^2
        cw = exp(-w)
        print(i, "     ")
        psi = psi + cw*BuildVBS(pathmats[i,:,:])
    end
    normalize!(psi)

    function Location(k::Int64,N::Int64)
        x = mod((k-1),N)+1
        y = ((k-1)÷N)+1
        return [x,y]
    end

    for i in 1:9, j in 1:9
        dvec = Location(i,3)-Location(j,3)
        push!(dvecnew, norm(dvec))
        cor = dot(psi,Bilin(i,j,9)*psi)
        push!(corfunc,cor)
        push!(abscorfunc,abs(cor))
    end
end

plot(dvecnew,abscorfunc,seriestype = :scatter, markeralpha = 0.15, ylabel = L"\left|\vec{S}_{i}\cdot\vec{S}_{j}\right|",xlabel = L"{\rm Separation} \  \left|\vec{x}_{i}-\vec{x}_{j}\right|",label = L"{\rm Exact} \ \theta" ,ytickfontsize = 10,xtickfontsize = 10,xguidefontsize = 15,yguidefontsize =15, fontfamily = "Times",legend = false, legendfontsize = 14)

savefig("correlationfunc2.pdf")



savearray = [dvecnew,corfunc,abscorfunc]
filename = string("alpha=0.1-beta=1.2averagedover2trialswith2x10^5.txt")
touch(filename)
open(filename, "w") do io
    writedlm(filename, savearray)
end










L=  9 #N^2
samples = 10^4
paths = walks(L,samples)
pathmats = WalkMats(paths)
avglength = 0
for i in 1:samples
    avglength = avglength+sum(Distances(paths[i],3))*(1/samples)*(1/L)
end

avglength
α = 0.1
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

function Location(k::Int64,N::Int64)
    x = mod((k-1),N)+1
    y = ((k-1)÷N)+1
    return [x,y]
end

for i in 1:9, j in 1:9
    dvec = Location(i,3)-Location(j,3)
    push!(dvecnew, norm(dvec))
    cor = dot(psi,Bilin(i,j,9)*psi)
    push!(corfunc,cor)
    push!(abscorfunc,abs(cor))
end






for i in 1:9
    print(dot(psi,Bilin(5,i,9) *psi),"   ")
end
dot(psi,Bilin(1,5,9) *psi)
Bilin(1,5,9)









L=  16 #N^2
samples = 500
paths = walks(L,samples)
pathmats = WalkMats(paths)
avglength = 0
for i in 1:samples
    avglength = avglength+sum(Distances(paths[i],4))*(1/samples)*(1/L)
end


avglength
α = 0.5
β = 1.9
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
cor = []
for i in 1:16
    correlation  = dot(psi,Bilin(9,i,16) *psi)
    push!(cor,correlation)
    print(correlation, "          ")
end
cor


test = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
Locations_xy([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16],4)
Distances(test,4)
