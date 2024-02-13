
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


    ## calculate mean field distribution (probability is NaN at the MF critical point)
    beta = 1.0 ./ ts
    heff = J*q*L0 .+ h
    pUp = exp.(beta .* heff) ./ (exp.(beta .* heff) .+ exp.(-beta .* heff))
    

    ## exclude critical temperature and take a simple sample
    idx = .!isnan.(pUp)
    pUp = pUp[idx]
    ts = ts[idx]
    

    ## Calculate the MF partition function
    logZm = -(0.5 * Nspins * q * J) * (L0.^2 ./ ts) .+ Nspins * (log.(cosh.(beta .* heff)) .+ log(2) )
    f1 = -(1 / Nspins) * logZm .* ts # feee energy mean field   
    f2 = -(1 / Nspins) * logZ .* ts # free energy trg
    # MakePlot(ts, f1, ts, f2, Tc, "TRG", "MF", "Free energies per unit spin", "T", "free energy per spin", missing)
    
    
    # Calculate constant for Accept-Reject algorithm
    htot = h .- heff
    
    C = 2*cosh.(beta .* heff)
    A = beta .* (2 * J .+ htot)

    logM = Nspins * (log.(C) .- (logZ ./ Nspins) .+ A)
    # logM = abs.(logM)

    display(scatter(ts, logM))

    ## sampling
    # nsamples = Int(ceil(exp(20)))
    # nsamples = 100
    # k = 30

    # @info ts[k]

    # samples = []
    # cnt = 0

    # while cnt < nsamples
    #     spins = sampleMF(Nspins, pUp[k])
    #     u = rand(1)[1]
    #     logP = -beta[k] * IsingEnergy(spins, J=J, h=h) - logZ[k]
    #     logPm = beta[k] * heff[k] * sum(spins) - Nspins * log.(C[k])
        
    #     log(u) < logP - logPm - logM[k] ? push!(samples,1) : push!(samples,0)

        # if log(u) < logP - logPm - logM[k] 
        #     push!(samples, spins)
        #     cnt += 1
        # end
   
    #     cnt += 1
    # end

    # sum(samples)

    
end