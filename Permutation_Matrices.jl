using LinearAlgebra
using Random
using SparseArrays

#=
N=5
L=N^2
MatInd = Vector{Int64}(1:L)
samples = 10000
Perms=[]
for i in 1:samples
    push!(Perms, shuffle(MatInd))
end
=#
#=
function walks(L::Int64,samples::Int64)
    MatInd = Vector{Int64}(1:L)
    Perms=[]
    for i in 1:samples
        push!(Perms, shuffle(MatInd))
    end
    return Perms
end
=#
function walks(L::Int64,samples::Int64)
    MatInd = Vector{Int64}(1:L)
    Perms=[]
    while length(Perms) < samples
        permutation = shuffle(MatInd)
        addval = true
        for i in 1:L
            if permutation[i] == i
                addval = false
            end
        end
        if addval
            push!(Perms, permutation)
        end
    end
    return Perms
end




function WalkMats(perms::Vector{Any})
    N = Int64(length(perms[1]))
    M = length(perms)
    walkmatrices = zeros(Int64,M,N,N)
    for i in 1:M
        for j in 1:N
            walkmatrices[i,j,perms[i][j]] = 1
        end
    end
    return walkmatrices
end

function Locations_xy(perm::Vector{Int64},N::Int64#=this is an NxN square lattice=#)
    step_dist = [] #where distances are stored for the individual steps
    xymap = []
    for i in 1: N^2
        x = (mod((perm[i]-1),N))+1
        y = ((perm[i]-1)÷N)+1
        push!(xymap, [x,y])
    end
    return xymap
end


function Distances(perm::Vector{Int64},N::Int64#=this is an NxN square lattice=#)
    positions = Locations_xy(perm,N)
    L = N^2
    distarray = []
    push!(distarray,norm(positions[1]-positions[L]))
    for i in 2:L
        push!(distarray,norm(positions[i]-positions[i-1]))
    end
    return Vector{Float64}(distarray)
end
#=
sum(Distances(walks(9,10)[1],3))
Locations_xy(walks(9,10)[3],3)
=#

#=
function Distances(positions,N::Int64)
    L = N^2
    distarray = []
    push!(distarray,norm(positions[1]-positions[L]))
    for i in 2:L
        push!(distarray,norm(positions[i]-positions[i-1]))
    end
    return Vector{Float64}(distarray)
end
=#


#=
distance = []
for i in 1:1000
    push!(distance,sum(Distances(walks(9,1000)[i],3)))
end
distance
min(distance...)

map = Locations_xy(walk,3)
=#
