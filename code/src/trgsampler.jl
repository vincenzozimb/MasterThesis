
## ---------- ACCEPT-REJECT SAMPLER FOR SAMPLING FROM THE TRG GIBBS DISTRIBUTION ---------- ##

using JLD

include("func.jl")
include("meanfield.jl")

let 

    ## Load data
    par = load("data/param.jld")
    res = load("data/resultsTrg.jld")
    

    J = get(par, "J", NaN)
    h = get(par, "h", NaN)

    Niter = get(par, "Niter", NaN)
    Nspins = 2^(Niter+1)
    @assert isinteger(sqrt(Nspins))
    
    ts = get(res, "ts", NaN)
    logZ = get(res, "logZ", NaN)
    

    ## Calculate MF magnetization solving the self consistent equation
    L0 = MFmagnetization(ts, J, h)
    

    ## Calculate the MF partition function
    heff = (J*q*L0 .+ h) ./ ts
    logZm = -(0.5 * Nspins * q * J) * (L0.^2 ./ ts) .+ Nspins * (log.(cosh.(heff)) .+ log(2) )

    f1 = -(1 / Nspins) * logZ .* ts # free energy trg
    f2 = -(1 / Nspins) * logZm .* ts # feee energy mean field   
    
    # MakePlot(ts, f1, ts, f2, Tc, "MF", "TRG", "Free energies per unit spin", "T", "free energy per spin", missing)
    
    # Calculate constant for Accept-Reject algorithm
    beta = 1 ./ ts
    
    logM = logZm .- logZ
    logM .+= (0.5 * J * Nspins * q) * beta .* (L0.^2)
    logM .+= (2 * Nspins * J) * beta
    
    htot = h*beta .- heff
    
    # logM .+= Nspins * sign.(htot) .* htot
       
    M = exp.(logM)

    scatter(ts, logM)

    
end