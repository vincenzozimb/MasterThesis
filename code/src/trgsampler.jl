
## ---------- ACCEPT-REJECT SAMPLER FOR SAMPLING FROM THE TRG GIBBS DISTRIBUTION ---------- ##

using JLD

include("func.jl")
include("meanfield.jl")


function AccRej(nsamples::Int, Nspins::Int, logM::Real, pMF::Real, logZ::Real, beta::Real, L0::Real, J=1.0, h=0.0)
    
    samples = []
    cnt = 0

    while cnt < nsamples
        spins = sampleMF_T(Nspins, pMF)
        u = rand(1)
        logP = -beta * IsingEnergy_T(spins, J=J, h=h) - logZ
        logPm = -beta * IsingEnergy_T(spins, J=0.0, h=(J*q*L0 + h))
        if log(u) < logP - logM - logPm
            push!(samples, spin)
            cnt += 1
        end
    end
    return samples
end


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
    beta = 1.0 ./ ts
    heff = (J*q*L0 .+ h) .* beta
    pUp = exp.(heff) ./ (exp.(heff) .+ exp.(-heff))

    logZm = -(0.5 * Nspins * q * J) * (L0.^2 ./ ts) .+ Nspins * (log.(cosh.(heff)) .+ log(2) )
    f1 = -(1 / Nspins) * logZm .* ts # feee energy mean field   
    f2 = -(1 / Nspins) * logZ .* ts # free energy trg
    # MakePlot(ts, f1, ts, f2, Tc, "TRG", "MF", "Free energies per unit spin", "T", "free energy per spin", missing)
    
    
    ## Calculate constant for Accept-Reject algorithm
    htot = h*beta .- heff
    
    C = 2*cosh.(heff)
    A = 2 * J .* beta .+ htot

    logM = Nspins * (log.(C) .- (logZ ./ Nspins) .+ A)
    

    nsamples = 10
    samples = AccRej(nsamples, Nspins, logM[10], pUp[10], logZ[10], beta[10], L0[10], J, h)

    # spin = sampleMF(Nspins, pUp)
    # E = IsingEnergy(spin, J=J, h=h) ./ Nspins
    # M = vec(mean(magn_sample, dims=1)) ./ Nspins
    # dM = vec(std(magn_sample, dims=1)) ./ Nspins
    # E = vec(mean(energy_sample, dims=1)) ./ Nspins
    # VarE = vec(var(energy_sample, dims=1)) ./ Nspins^2

    # MakePlot(ts, E, NaN, NaN, Tc, "MF", "", "Energy", "T", "E", missing)
    # MakeErrorPlot(ts, E, sqrt.(VarE), Tc, "MF", "Energy", "T", "E", missing)

    
end