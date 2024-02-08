## ---------- Run TRG algorithm ---------- ##

using JLD

include("trg.jl")
include("func.jl")

let 

    ## Physical parameters
    
    # bond dimension (max 15 on my laptop)
    Dcut = 6
    # number of iterations
    Niter = 15
    # coupling constant
    J = 1.0
    # external field
    h = 0.01


    ## Simulation parameters

    # temperature range
    ts = 0.1:0.15:10.0
    # βJ range
    ks = J ./ ts
    

    ## TRG run
    welcome = "=====TRG======" 
    @show welcome
    logZ = []
    logZh = []
    cnt = size(ks,1)
    for k in ks
        A = tensor_chess(k, 0.0)
        y = trg(A, Dcut, Niter)
        push!(logZ,y)
        A = tensor_chess(k, h)
        y = trg(A, Dcut, Niter)
        push!(logZh,y)
        print("\r count=$cnt")
        flush(stdout)
        cnt -= 1
    end
    println()
    
    
    ## Save results
    # save("results.jld")


    ## Calculate results
    
    # free energies
    F = - ts .* logZ / 2^(Niter+1);
    Fh = -ts .* logZh / 2^(Niter+1);
    Fexact = ising_free_energy.(1.0 ./ ts, J)
    
    # relative error in free energy
    re = abs.(F - Fexact) ./ abs.(Fexact) 

    #  specific heat and magnetization
    M = -(Fh - F) ./ (h * ts)
    C = SpecificHeat(F,ts)
    
    tp = 0.1:0.01:10.0
    Fexactp = ising_free_energy.(1.0 ./ tp, J)
    Cexact = SpecificHeat(Fexactp,tp)
    Mexact = ising_magnetization.(1.0 ./ tp)


    ## Plots
    MakePlot(ts, F, ts, Fexact, "TRG", "Exact", "Free Energy per site", "T", "F", "FreeEnergy.png")
    MakePlot(ts, re, NaN, NaN, "TRG", "", "Free energy relative error", "T", "ϵ", "RelativeError.png")
    MakePlot(ts, C, tp, Cexact, "TRG", "Exact", "Specific heat", "T", "C", "SpecificHeat.png")
    MakePlot(ts, M, tp, Mexact, "TRG", "Exact", "Magnetization", "T", "M", "Magnetization.png")
    
end