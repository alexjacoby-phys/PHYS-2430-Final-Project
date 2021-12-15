using LinearAlgebra
using SparseArrays

projup = sparse([1 0 0 ; 0 0 0 ; 0 0 0])
projzero = sparse([0 0 0 ; 0 1 0 ; 0 0 0])
projdown = sparse([0 0 0 ; 0 0 0 ; 0 0 1])
PJ = [projup,projzero,projdown]
trup = sparse([1 0 0])
trzero = sparse([0 1 0])
trdown = sparse([0 0 1])
TR = [trup,trzero,trdown]
Id3 = sparse(I,3,3)


function PJOP(k::Int64#=pick which state, 1 = up, 2 = zero, 3 = down=#,n::Int64,L::Int64)
    projop = PJ[k]
    oparray = []
    for i in 1:L
        if i != n
            push!(oparray,Id3)
        elseif i == n
            push!(oparray,projop)
        else
            return "error"
        end
    end
    return(kron(oparray...))
end
function TROP(k::Int64#=pick which state, 1 = up, 2 = zero, 3 = down=#,n::Int64,L::Int64)
    tracop = TR[k]
    oparray = []
    for i in 1:L
        if i != n
            push!(oparray,Id3)
        elseif i == n
            push!(oparray,tracop)
        else
            return "error"
        end
    end
    return kron(oparray...)
end


#=
psi = sparse([1 1 1])
adjoint(psi)
kron(sparse(psi),sparse(adjoint(psi)))
=#



function Collapse(k::Int64 #=which state projected upon=#,n::Int64, L::Int64,psi)
    psi = PJOP(k,n,L)*psi
    return psi
end

function DM_Build(psi)
    return kron(psi,adjoint(psi))
end


function Sparse_DM_Build(psi)
    return kron(sparse(psi),adjoint(sparse(psi)))
end

function Site_Trace(rho,k::Int64,L::Int64)
    retrho = TROP(1,k,L)*rho*adjoint(TROP(1,k,L))+TROP(2,k,L)*rho*adjoint(TROP(2,k,L))+TROP(3,k,L)*rho*adjoint(TROP(3,k,L))
    return retrho
end

function Partial_Trace(rho,sites0::Vector{Int64},L::Int64)
    sites = sort(sites0, rev = true)
    N = L
    k = length(sites)
    DM = rho
    for i in sites
        DM = Site_Trace(DM,i,N)
        N= N-1
    end
    return DM
end

function SvN(svalues::Vector{Float64},cutoff::Float64)
    schmidtvec = real(droptol!(sparse(svalues),cutoff))
    S=0
    for i in 1:length(schmidtvec)
        S = S+ schmidtvec[i]*log(schmidtvec[i])
    end
    return S
end

#function Sparse_TR_DM(sites::Vector{Int64}) this is a test
