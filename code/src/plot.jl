using JLD

include("func.jl")

let 

    ## load result from TRG algorithm
    d = load("data/results.jld")

    J = get(d, "J", NaN)
    h = get(d, "h", NaN)
    ts = get(d, "ts", NaN)
    F = get(d, "F", NaN)
    Fh = get(d, "Fh", NaN)

    
    ## Calculate observables
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