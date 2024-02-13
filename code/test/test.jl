using Plots
using Statistics

# using StatsPlots

include("../src/func.jl")

let 

    u = rand(10000)
    histogram(u)    

end