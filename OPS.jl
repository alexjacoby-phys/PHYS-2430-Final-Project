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



function Biquad(n::Int64,m::Int64,L::Int64)
    opsarray11 = []
    opsarray21 = []
    opsarray31 = []
    opsarray12 = []
    opsarray22 = []
    opsarray32 = []
    opsarray13 = []
    opsarray23 = []
    opsarray33 = []
    if m == n
        for i in 1:L
            push!(opsarray11,id)
            push!(opsarray21,id)
            push!(opsarray31,id)
            push!(opsarray12,id)
            push!(opsarray22,id)
            push!(opsarray32,id)
            push!(opsarray13,id)
            push!(opsarray23,id)
            push!(opsarray33,id)
        end
        opsarray11[n] = (sz^2)*(sz^2)
        opsarray21[n] = (sy^2)*(sz^2)
        opsarray31[n] = (sx^2)*(sz^2)
        opsarray12[n] = (sz^2)*(sy^2)
        opsarray22[n] = (sy^2)*(sy^2)
        opsarray32[n] = (sx^2)*(sy^2)
        opsarray13[n] = (sz^2)*(sx^2)
        opsarray23[n] = (sy^2)*(sx^2)
        opsarray33[n] = (sx^2)*(sx^2)
        return kron(opsarray11...)+kron(opsarray21...)+kron(opsarray31...)+kron(opsarray12...)+kron(opsarray22...)+kron(opsarray32...)+kron(opsarray13...)+kron(opsarray23...)+kron(opsarray33...)
    else
        for i in 1:L
            if (i != n) & (i != m)
                push!(opsarray11,id)
                push!(opsarray21,id)
                push!(opsarray31,id)
                push!(opsarray12,id)
                push!(opsarray22,id)
                push!(opsarray32,id)
                push!(opsarray13,id)
                push!(opsarray23,id)
                push!(opsarray33,id)
            elseif i == n # these entries are wrong-- these should be bilinear in these operators
                push!(opsarray11,sz^2)
                push!(opsarray21,sz*sy)
                push!(opsarray31,sz*sx)
                push!(opsarray12,sy*sz)
                push!(opsarray22,sy^2)
                push!(opsarray32,sy*sx)
                push!(opsarray13,sx*sz)
                push!(opsarray23,sx*sy)
                push!(opsarray33,sx^2)
            elseif i == m
                push!(opsarray11,sz^2)
                push!(opsarray21,sz*sy)
                push!(opsarray31,sz*sx)
                push!(opsarray12,sy*sz)
                push!(opsarray22,sy^2)
                push!(opsarray32,sy*sx)
                push!(opsarray13,sx*sz)
                push!(opsarray23,sx*sy)
                push!(opsarray33,sx^2)
            end
        end
        return kron(opsarray11...)+kron(opsarray21...)+kron(opsarray31...)+kron(opsarray12...)+kron(opsarray22...)+kron(opsarray32...)+kron(opsarray13...)+kron(opsarray23...)+kron(opsarray33...)
    end
end
