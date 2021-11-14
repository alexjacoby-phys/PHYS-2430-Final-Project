using LinearAlgebra
using Arpack
using SparseArrays

#=\

These are no longer accurate!
sites are linked in the following pattern in the 3x3 case
. _ . _  .
|. _ . _ .
. _ . _ .|
where these flow from the bottom left to top right 1-9

similarly in the 4x4 case:
. _ . _ . _  .
. _ . _ . _ .|
|. _ . _ . _ .
. _ . _ . _ .|
with the same numbering scheme

here the numbering scheme denotes the ordering in the Kronecker product array

=#

id = (1+0*im)*sparse([1 0 0; 0 1 0; 0 0 1])
sz = (1+0*im)*sparse([1 0 0 ; 0 0 0; 0 0 -1])
sx = (1+0*im)*sparse((1/sqrt(2))*[0 1 0; 1 0 1; 0 1 0 ])
sy = (1+0*im)*sparse((1/sqrt(2))*[0 -im 0; im 0 -im ; 0 im 0])
up = (1+0*im)*sparse([1,0,0])
zero = (1+0*im)*sparse([0,1,0])
down = (1+0*im)*sparse([0,0,1])
AKLT = zeros(Complex{Float64},2,3,2)
AKLT[1,:,1] = (1/sqrt(2))*zero
AKLT[2,:,2] = -(1/sqrt(2))*zero
AKLT[1,:,2] = -up
AKLT[2,:,1] = down


#= test = zeros(Int64,2,4)
test[1,:] = [1,1,1,1]
test[2,:] = [2,2,2,2]
test
run with test to recover the wavefunction ψ = [Λ_{1}]^{1}_{2} ⊗ [Λ_{2}]^{1}_{2} ⊗ [Λ_{3}]^{1}_{2} ⊗ [Λ_{4}]^{1}_{2}    =#

function statbuild(statdex::Array{Int64,2})
    retarray = []
    for i in 1:length(statdex[1,:])
        push!(retarray,sparse(AKLT[statdex[1,i],:,statdex[2,i]]))
    end
    return kron(retarray...)
end
#=
function statbuild(statdex::Array{Int64,2})
    retarray = []
    for i in 1:length(statdex[1,:])
        push!(retarray,AKLT[statdex[1,i],:,statdex[2,i]])
    end
    return sparse(kron(retarray...))
end
=#
#=

VBS wavefunctions are built with a matrix that links up indices since we must contract them non-locally in the 2d to 1d mapping. Here we will consider the row index of the MPS to be the row of the contraction matrix (denoted ContraMat, which specifies which indices are to tie with which). A contraction matrix indicates a specific tiling of singlet bonds on the 2d lattice.

For example, if we consider the singlet pattern of the upper index of site 1 contracting with the lower index of site 3, the upper index of site 3 contracting with the lower index of site 2, and the upper index of site 2 contracting with the lower index of site 1 this will look like
[0 0 1; 1 0 0; 0 1 0]

=#


function Contractions(ContraMat::Array{Int64})
    # this bit finds which indices should be contracted. Row indices of MPS are placed in always the first of the dyads
    nz = findnz(sparse(ContraMat))
    IndTies = []
    for i in 1:length(nz[1])
        push!(IndTies,[nz[1][i],nz[2][i]])
    end
    return IndTies
    #=this returns a list of 2 entry vectors where the first entry is the upper index site and the second entry is the lower index site that are to be contracted against one another=#
end

#= For a check, run Contractions([0 0 1; 1 0 0; 0 1 0]) to see that we recover the example that I mention above. It should produce the list Vector{Any} with 3 elements, [[2,1], [3,2], [1,3]] =#

#=Now we have to do the dirty work of inputting such a list of contractions to get back a list of indices to input into statbuild for each VBS configuration specified by a contramat the input of this function is just indties, the output of contractions (labeled CTRS)=#

#CTRS = Contractions([0 0 1; 1 0 0; 0 1 0])

function IndexBuilder(CTRS::Vector{Any})
#=there is some dumb bug in doing this with nested arrays so you just have to define it like this-- the first index chooses the upper (1) vs lower (2) spin index on on a given site, the second index is the site (1-L), and the third index is which index value configuration within the VBS we are in  (goes 1-2^L for all combinations of indices). The actual value is the spin index.=#
    L=length(CTRS)
    indices_for_buildstat = zeros(Int64,2,L,2^L)
    for i in 1:length(CTRS), j in 1:2^(length(CTRS))
        #print(mod((j-1)÷2^(length(CTRS)-i),2))
        ContrSite = CTRS[i]
        if mod((j-1)÷2^(length(CTRS)-i),2) == 0
            indices_for_buildstat[1,ContrSite[1],j] = indices_for_buildstat[2,ContrSite[2],j] = 1
        elseif mod((j-1)÷2^(length(CTRS)-i),2) == 1
            indices_for_buildstat[1,ContrSite[1],j] = indices_for_buildstat[2,ContrSite[2],j] = 2
        else
            return "error"
        end
    end
    return indices_for_buildstat
end
#=
In the spirit of testable code, run the below thing to see this is working

CTRS = Contractions([0 0 1; 1 0 0; 0 1 0])

L=length(CTRS)
indices_for_buildstat = zeros(Int64,2,L,2^L)

for i in 1:length(CTRS), j in 1:2^(length(CTRS))
    print(mod((j-1)÷2^(length(CTRS)-i),2))
    ContrSite = CTRS[i]
    if mod((j-1)÷2^(length(CTRS)-i),2) == 0
        indices_for_buildstat[1,ContrSite[1],j] = indices_for_buildstat[2,ContrSite[2],j] = 1
    elseif mod((j-1)÷2^(length(CTRS)-i),2) == 1
        indices_for_buildstat[1,ContrSite[1],j] = indices_for_buildstat[2,ContrSite[2],j] = 2
    else
        return "error"
    end
end
indices_for_buildstat


it is also of non-trivial value to run a sort of conjoined test case with the previous function just to make sure everything is working smoothly together
=#

function BuildVBS(ContraMat::Array{Int64,2})
    Contr = Contractions(ContraMat)
    Indices = IndexBuilder(Contr)
    L = length(ContraMat[:,1])
    psi = zeros(Float64,3^L)
    for i in 1:2^L
        psi = psi + sparse(statbuild(Indices[:,:,i]))
    end
    return normalize(psi)
end



#ContraMat = [0 1 0 0 ; 0 0 0 1 ; 1 0 0 0 ; 0 0 1 0]
#BuildVBS(ContraMat)
