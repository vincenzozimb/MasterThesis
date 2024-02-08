## ---------- MEAN FIELD THEORY FOR SQUARE LATTICE ISING MODEL ---------- ##

using Plots
using JLD

include("func.jl")


## function for self consistent equation 
function mf(L::Real, β::Real; h::Real=0.0, J::Real=1.0, q::Int=4)
    return L - tanh(β * (h + q*J*L) )
end


## derivative of the function for self consistent equation
function mfder(L::Real, β::Real; h::Real=0.0, J::Real=1.0, q::Int=4)
    return 1.0 - β*J*q / (cosh( β * (h + q*J*L) )^2)
end


let 
    
    ## load data
    d = load("data/results.jld")
    
    J = get(d, "J", NaN)
    h = get(d, "h", NaN)
    ts = get(d, "ts", NaN)
    
    
    ## calculate zeros
    x0 = 1.0
    Tcmf = 4.0 # mean field critical temperature
    
    L0 = []
    for t in ts
        t != Tcmf ? y = NewtonRaphson(mf, x0, mfder, (1.0/t,)) : y = NaN 
        push!(L0,y)
    end


    ## calculate mean field distribution
    q = 4
    heff = (J*q*L0 .+ h) ./ ts
    pUp = exp.(heff) ./ (exp.(heff) .+ exp.(-heff))
    # probDn = exp.(-heff) ./ (exp.(heff) .+ exp.(-heff))
    # @info probUp .+ probDn

    # @show ts[26] ts[27] ts[28]
    @info "probabilities" pUp[26] pUp[27] pUp[28]

    
    nsamples = 20000
    samples = [[rand(Bernoulli(p)) for _ in 1:nsamples] for p in filter(!isnan, pUp)]
    # spins = 2 .* samples .- 1 


    # samples = [rand(Bernoulli(pUp[60])) for _ in 1:nsamples]

    # M = sum(spins) / nsamples
    # @info M

    # samples = [rand(Bernoulli(p)) for p in filter(!isnan, pUp)]
    # nothing

end