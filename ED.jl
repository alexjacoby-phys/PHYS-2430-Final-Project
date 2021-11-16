include("OPS.jl")
using LinearAlgebra
using SparseArrays
using Arpack

N=3
L=N^2
Ham1 = spzeros(Float64,3^L,3^L)
Ham2 = spzeros(Float64,3^L,3^L)


Ham1 = Bilin(1,2,9)+Bilin(2,3,9)+Bilin(4,5,9)+Bilin(5,6,9)+Bilin(7,8,9)+Bilin(8,9,9)+Bilin(1,4,9)+Bilin(2,5,9)+Bilin(3,6,9)+Bilin(4,7,9)+Bilin(5,8,9)+Bilin(6,9,9)

Ham2 = Biquad(1,2,9)+Biquad(2,3,9)+Biquad(4,5,9)+Biquad(5,6,9)+Biquad(7,8,9)+Biquad(8,9,9)+Biquad(1,4,9)+Biquad(2,5,9)+Biquad(3,6,9)+Biquad(4,7,9)+Biquad(5,8,9)+Biquad(6,9,9)

function H(θ::Float64)
    H = cos(θ)*Ham1+sin(θ)*Ham2
    return H
end




#ED= eigs(H, nev=2, which = :SR, maxiter = 100)

#psi0= ED[2][:,1]

#dot(psi,psi0)

#dot(psi,H*psi)
