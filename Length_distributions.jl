include("StateBuilder.jl")
include("Permutation_Matrices.jl")
include("InfoTools.jl")
include("Ops.jl")
using Plots



N=20
L=  N^2
samples = 10000
paths = walks(L,samples)
pathmats = WalkMats(paths)

lengths = Vector{Float64}([])
avglength = 0
for i in 1:samples
    dist = sum(Distances(paths[i],N))
    push!(lengths,dist)
    avglength = avglength + dist*(1/samples)
end
avglength
dat = (1/avglength)*sort(lengths)

datlength = length(dat)
spaces = [0.]
for i in 2:datlength
    push!(spaces,dat[i]-dat[i-1])
end

cleandat = Vector{Float64}([])
cleaninvsp = Vector{Float64}([])

for i in 1:datlength
    if !iszero(round(spaces[i],digits = 6))
        push!(cleandat,dat[i])
        push!(cleaninvsp,(spaces[i])^-1)
    end
end
cleaninvsp = cleaninvsp/max(cleaninvsp...)

plot(cleandat, cleaninvsp, legend = false,seriestype = :scatter, alpha = 0.02)
