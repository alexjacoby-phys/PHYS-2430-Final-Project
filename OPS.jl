using SparseArrays
using LinearAlgebra
id = (1+0*im)*sparse([1 0 0; 0 1 0; 0 0 1])
sz = (1+0*im)*sparse([1 0 0 ; 0 0 0; 0 0 -1])
sx = (1+0*im)*sparse((1/sqrt(2))*[0 1 0; 1 0 1; 0 1 0 ])
sy = (1+0*im)*sparse((1/sqrt(2))*[0 -im 0; im 0 -im ; 0 im 0])

function Bilin(n::Int64,m::Int64,L::Int64)
    opsarray1 = []
    opsarray2 = []
    opsarray3 = []
    if m == n
        for i in 1:L
            push!(opsarray1,id)
            push!(opsarray2,id)
            push!(opsarray3,id)
        end
        opsarray1[n] = sz*sz
        opsarray2[n] = sy*sy
        opsarray3[n] = sx*sx
        return kron(opsarray1...)+kron(opsarray2...)+kron(opsarray3...)
    else
        for i in 1:L
            if (i != n) & (i != m)
                push!(opsarray1,id)
            elseif i == n
                push!(opsarray1,sz)
            elseif i == m
                push!(opsarray1,sz)
            end

            if (i != n) & (i != m)
                push!(opsarray2,id)
            elseif i == n
                push!(opsarray2,sy)
            elseif i == m
                push!(opsarray2,sy)
            end

            if (i != n) & (i != m)
                push!(opsarray3,id)
            elseif i == n
                push!(opsarray3,sx)
            elseif i == m
                push!(opsarray3,sx)
            end
        end
        return kron(opsarray1...)+kron(opsarray2...)+kron(opsarray3...)
    end
end
