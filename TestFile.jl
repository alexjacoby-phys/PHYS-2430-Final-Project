include("StateBuilder.jl")
include("Permutation_Matrices.jl")
include("InfoTools.jl")
include("Ops.jl")
include("ED.jl")
using Plots


#first group of small tests.
mat = [0 1 ; 1 0]
CTRS = Contractions(mat)
IndexBuilder(CTRS)

mat = [0 0 1 ; 1 0 0 ; 0 1 0]
CTRS = Contractions(mat)
IndexBuilder(CTRS)

mat = [0 0 0 1 ; 1 0 0 0 ; 0 1 0 0 ; 0 0 1 0]
CTRS = Contractions(mat)
IndexBuilder(CTRS)


# second group of small tests
mat = [0 0 1 ; 1 0 0 ; 0 1 0]
psi = normalize(BuildVBS(mat))
Ham13 = (1/2)*Bilin(1,3,3)+(1/6)*Biquad(1,3,3)
Ham21 = (1/2)*Bilin(2,1,3)+(1/6)*Biquad(2,1,3)
Ham32 = (1/2)*Bilin(2,1,3)+(1/6)*Biquad(2,1,3)
dot(psi, Ham13*psi)
dot(psi, Ham21*psi)
dot(psi, Ham32*psi)

mat  =  [0 0 0 1 ; 1 0 0 0 ; 0 1 0 0 ; 0 0 1 0]
psi = normalize(BuildVBS(mat))
Ham14 = (1/2)*Bilin(1,4,4)+(1/6)*Biquad(1,4,4)
Ham21 = (1/2)*Bilin(2,1,4)+(1/6)*Biquad(2,1,4)
Ham32 = (1/2)*Bilin(3,2,4)+(1/6)*Biquad(3,2,4)
Ham43 = (1/2)*Bilin(4,3,4)+(1/6)*Biquad(4,3,4)
dot(psi, Ham14*psi)
dot(psi, Ham21*psi)
dot(psi, Ham32*psi)
dot(psi, Ham43*psi)



# third group of small tests
mat = [0 0 1 ; 1 0 0 ; 0 1 0]
psi = normalize(BuildVBS(mat))
DM = Sparse_DM_Build(psi)
eigvals(Matrix(Partial_Trace(DM, [1],3)))
eigvals(Matrix(Partial_Trace(DM, [2],3)))
eigvals(Matrix(Partial_Trace(DM, [3],3)))


mat = [0 1 0 ; 0 0 1 ; 1 0 0]
psi = normalize(BuildVBS(mat))
DM = Sparse_DM_Build(psi)
eigvals(Matrix(Partial_Trace(DM, [1],3)))
eigvals(Matrix(Partial_Trace(DM, [2],3)))
eigvals(Matrix(Partial_Trace(DM, [3],3)))

mat = [0 0 0 1 ; 1 0 0 0 ; 0 1 0 0 ; 0 0 1 0]
psi = normalize(BuildVBS(mat))
DM = Sparse_DM_Build(psi)
eigvals(Matrix(Partial_Trace(DM, [1],4)))
eigvals(Matrix(Partial_Trace(DM, [2],4)))
eigvals(Matrix(Partial_Trace(DM, [3],4)))

mat = [0 1 0 0 ; 0 0 0 1 ; 1 0 0 0 ; 0 0 1 0]
psi = normalize(BuildVBS(mat))
DM = Sparse_DM_Build(psi)
eigvals(Matrix(Partial_Trace(DM, [1],4)))
eigvals(Matrix(Partial_Trace(DM, [2],4)))
eigvals(Matrix(Partial_Trace(DM, [3],4)))






#fourth group of small tests

paths = [2,4,1,3]
Distances(paths,2)
sum(Distances(paths,2))



paths = [4,1,2,3]
Distances(paths,2)
sum(Distances(paths,2))


paths = [3,4,5,1,7,8,9,2,6]
Distances(paths,3)
sum(Distances(paths,3))



#first at scale check

L=  4 #N^2
samples = 10^5
paths = walks(L,samples)
pathmats = WalkMats(paths)
avglength = 0
for i in 1:samples
    avglength = avglength+sum(Distances(paths[i],2))*(1/samples)*(1/L)
end

avglength
psi = zeros(Complex{Float64},3^L)
for i in 1:samples
    dist = sum(Distances(paths[i],2))/L
    psi = psi + BuildVBS(pathmats[i,:,:])
    #print("x")
end
normalize!(psi)

DM = Sparse_DM_Build(psi)
eigvals(Matrix(Partial_Trace(DM, [1],4)))
eigvals(Matrix(Partial_Trace(DM, [2],4)))
eigvals(Matrix(Partial_Trace(DM, [3],4)))
eigvals(Matrix(Partial_Trace(DM, [4],4)))
eigvals(Matrix(Partial_Trace(DM, [2,3,4],4)))
eigvals(Matrix(Partial_Trace(DM, [1,3,4],4)))
eigvals(Matrix(Partial_Trace(DM, [1,2,4],4)))
eigvals(Matrix(Partial_Trace(DM, [2,3,4],4)))


#second at scale check

L=  4 #N^2
samples = 10^5
paths = walks(L,samples)
pathmats = WalkMats(paths)
avglength = 0
for i in 1:samples
    avglength = avglength+sum(Distances(paths[i],2))*(1/samples)*(1/L)
end

avglength
α = 0.1
β = 1.05
psi = zeros(Complex{Float64},3^L)
for i in 1:samples
    dist = sum(Distances(paths[i],2))/L
    w = (0.5*(dist-β)/α)^2
    cw = exp(-w)
    print(cw)
    psi = psi + cw*BuildVBS(pathmats[i,:,:])
end
normalize!(psi)

DM = Sparse_DM_Build(psi)
eigvals(Matrix(Partial_Trace(DM, [1],4)))
eigvals(Matrix(Partial_Trace(DM, [2],4)))
eigvals(Matrix(Partial_Trace(DM, [3],4)))
eigvals(Matrix(Partial_Trace(DM, [4],4)))
eigvals(Matrix(Partial_Trace(DM, [2,3,4],4)))
eigvals(Matrix(Partial_Trace(DM, [1,3,4],4)))
eigvals(Matrix(Partial_Trace(DM, [1,2,4],4)))
eigvals(Matrix(Partial_Trace(DM, [2,3,4],4)))

eigvals(Matrix(Partial_Trace(DM, [1,2],4)))
eigvals(Matrix(Partial_Trace(DM, [1,3],4)))




L=  4 #N^2
samples = 10^6
paths = walks(L,samples)
pathmats = WalkMats(paths)
avglength = 0
for i in 1:samples
    avglength = avglength+sum(Distances(paths[i],2))*(1/samples)*(1/L)
end

avglength
α = 0.1
β = 1.05
psi = zeros(Complex{Float64},3^L)
for i in 1:samples
    dist = sum(Distances(paths[i],2))/L
    w = (0.5*(dist-β)/α)^2
    cw = exp(-w)
    print(cw)
    psi = psi + cw*BuildVBS(pathmats[i,:,:])
end
normalize!(psi)

DM = Sparse_DM_Build(psi)
eigvals(Matrix(Partial_Trace(DM, [1,2],4)))
eigvals(Matrix(Partial_Trace(DM, [1,3],4)))



#third at scale check
L=  9 #N^2
samples = 10^4
paths = walks(L,samples)
pathmats = WalkMats(paths)
psi = zeros(Complex{Float64},3^L)
for i in 1:samples
    #print(cw)
    psi = psi + BuildVBS(pathmats[i,:,:])
end
normalize!(psi)

DM = Sparse_DM_Build(psi)
eigvals(Matrix(Partial_Trace(DM, [1],9)))
eigvals(Matrix(Partial_Trace(DM, [2],9)))
eigvals(Matrix(Partial_Trace(DM, [3],9)))
eigvals(Matrix(Partial_Trace(DM, [4],9)))
eigvals(Matrix(Partial_Trace(DM, [2,3,4,5,6,7,8,9],9)))
eigvals(Matrix(Partial_Trace(DM, [1,2],9)))
eigvals(Matrix(Partial_Trace(DM, [1,3],9)))
eigvals(Matrix(Partial_Trace(DM, [1,4],9)))




#fifth at scale check
L=  9 #N^2
samples = 10^6
paths = walks(L,samples)
pathmats = WalkMats(paths)
avglength = 0
for i in 1:samples
    avglength = avglength+sum(Distances(paths[i],3))*(1/samples)*(1/L)
end

avglength
α = 0.1
β = 1.4
psi = zeros(Complex{Float64},3^L)
for i in 1:samples
    dist = sum(Distances(paths[i],3))/L
    w = (0.5*(dist-β)/α)^2
    cw = exp(-w)
    #print(cw)
    psi = psi + cw*BuildVBS(pathmats[i,:,:])
end
normalize!(psi)

DM = Sparse_DM_Build(psi)
eigvals(Matrix(Partial_Trace(DM, [1],9)))
eigvals(Matrix(Partial_Trace(DM, [2],9)))
eigvals(Matrix(Partial_Trace(DM, [3],9)))
eigvals(Matrix(Partial_Trace(DM, [4],9)))
eigvals(Matrix(Partial_Trace(DM, [2,3,4,5,6,7,8,9],9)))
eigvals(Matrix(Partial_Trace(DM, [1,2],9)))
eigvals(Matrix(Partial_Trace(DM, [1,3],9)))
eigvals(Matrix(Partial_Trace(DM, [1,4],9)))
eigvals(Matrix(Partial_Trace(DM, [1,7],9)))




### everything after this is scratch.


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



DM = Sparse_DM_Build(psi)
eigvals(Matrix(Partial_Trace(DM, [1,2,3,4,5,6],9)))


for i in 1:9
    print(dot(psi,Bilin(5,i,9) *psi),"   ")
end
dot(psi,Bilin(1,5,9) *psi)
Bilin(1,5,9)

DM = Sparse_DM_Build(normalize(Collapse(1,5,9,psi)))
Partial_Trace(DM,[1,2,3,4,6,7,8,9],9)
#=

Distances(paths[1],3)
paths[1]


Locations_xy(paths[1],3)
=#


#=
aklt0 = normalize((kron(AKLT[1,:,1],AKLT[1,:,1],)+kron(AKLT[1,:,2],AKLT[2,:,1])))
adjoint(aklt0)*AKLTHAM*aklt0
AKLTHAM = (1/6)*Bilin(1,2,2)*Bilin(1,2,2)+(1/2)*Bilin(1,2,2)+(1/3)*kron(Id3,Id3)
always good to check
=#
#=
DM = Sparse_DM_Build(normalize!(Collapse(1,5,9,psi)))
for i in 0:7
    DM = Site_Trace(DM,9-i,9-i)
end

DM
=#
