
## ---------- TRG GIBBS DISTRIBUTION vs MEAN FIELD ---------- ##

using JLD
# using Statistics

include("func.jl")


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
    logZh = get(res, "logZh", NaN)
    
    
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
    
    
    # Calculate constant for Accept-Reject algorithm
    htot = h .- heff
    
    C = 2*cosh.(beta .* heff)
    A = beta .* (2 * J .+ abs.(htot))
    
    logM = Nspins * (log.(C) .- (logZ ./ Nspins) .+ A)
    
    # create trg images subfolder
    images_path = pwd() * "/images/stat/"
    if !isdir(images_path)
        mkdir(images_path)
    end
    
    
    # mean field vs TRG free energy
    Make2Plot(ts, f1, ts, f2, Tc, "TRG", "MF", "Free energies per unit spin", L"T", L"F/N", images_path, "TRGvsMF.png")
    

    # log-acceptance probability
    plt = scatter(ts, logM, ms=2, label=missing, yaxis=:log)
    vline!([Tc], line=:red, label=L"T_c")
    vline!([4], line=:green, label=L"T_c^{MF}")
    title!("Log-Acceptance Rate Gibbs vs MF")
    xlabel!(L"T")
    ylabel!(L"\log (M)")
    savefig(joinpath(images_path, "LogAcceptanceRate.png"))
    

    ## KL divergence of Gibbs from Mean Field

    # sample from mean field
    nsamples = 10


    SampleAver1 = zeros(Real, length(ts), nsamples)
    SampleAver2 = zeros(Real, length(ts), nsamples)
    for n in 1:nsamples
        samples = sampleMF.(Nspins,pUp)
        for t in 1:length(ts)
            SampleAver1[t,n] = IsingEnergy(samples[t], J=J, h=0.0)
            SampleAver2[t,n] = sum(samples[t])
        end
    end
    
    Aver1 = vec(mean(SampleAver1, dims=2))
    Aver2 = vec(mean(SampleAver2, dims=2))
    dAver1 = vec(std(SampleAver1, dims=2))
    dAver2 = vec(std(SampleAver2, dims=2))


    Dm = -beta .* (J * Aver1 .+ htot .* Aver2) .- (Nspins * log.(C)) .+ logZ
    dDm = beta .* sqrt.((J^2 * dAver1.^2) .+ (htot.^2 .* dAver2.^2))
    
    # MakeErrorPlot(ts, Dm, 100*dDm, Tc, "KL div", "KL divergence from MF", "T", "D_{KL}", missing, missing)
    Make2Plot(ts, Dm, NaN, NaN, [Tc, 4], "KL div", "", "KL divergence from MF", L"T", L"D_{KL}", images_path, "KLm.png")
    MakePlotLog(ts, Dm, [Tc, 4], "", "KL divergence from MF", L"T", L"D_{KL}", images_path, "KLmLog.png")


    ## KL divergence of Mean Field from Gibbs

    M = (logZh .- logZ) ./ (h * beta)
    Aver = diff(logZ) ./ diff(beta)
    Aver ./= J
    push!(Aver, NaN)

    Dm = -beta .* (-J * Aver .- htot .* M) .+ (Nspins * log.(C)) .- logZ

    scatter(ts, Dm)
    Make2Plot(ts, Dm, NaN, NaN, [Tc, 4], "KL div", "", "KL divergence from Gibbs", L"T", L"D_{KL}", images_path, "KL.png")
    MakePlotLog(ts, Dm, [Tc, 4], "KL div", "KL divergence from Gibbs", L"T", L"D_{KL}", images_path, "KLlog.png")

    # prova = []
    # for k in 1:length(ts)
    #     push!(prova, IsingEnergy(sampleMF(Nspins,pUp[k]), J=J, h=htot[k]))
    # end
    # scatter(ts, prova)


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