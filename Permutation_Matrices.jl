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

function walks(L::Int64,samples::Int64)
    MatInd = Vector{Int64}(1:L)
    Perms=[]
    for i in 1:samples
        push!(Perms, shuffle(MatInd))
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
