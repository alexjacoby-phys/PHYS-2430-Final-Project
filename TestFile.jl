include("StateBuilder.jl")
include("Permutation_Matrices.jl")
include("InfoTools.jl")
include("Ops.jl")





L=  9 #N^2
samples = 1000

psi = zeros(Complex{Float64},3^L)
paths = walks(L,samples)
pathmats = WalkMats(paths)

avglength = 0
for i in 1:samples
    avglength = avglength+sum(Distances(paths[i],3))*(1/samples)*(1/L)
end
avglength
α = 0.2
β = 2
for i in 1:samples
    dist = sum(Distances(paths[i],3))/L
    w = ((dist-β)/α)^2
    cw = ((avglength-β)/α)^2
    A = exp(cw-w)
    #print(A,"     ")
    psi = psi + A*BuildVBS(pathmats[i,:,:])
    print("x")
end
normalize!(psi)

DM = Sparse_DM_Build(psi)
eigen(Matrix(Partial_Trace(DM, [4,5,6,7,8,9],9)))


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
