using JLD

include("func.jl")

let 

    ## load result from TRG algorithm
    par = load("data/param.jld")
    res = load("data/resultsTrg.jld")

    J = get(par, "J", NaN)
    h = get(par, "h", NaN)
    ts = get(res, "ts", NaN)
    F = get(res, "F", NaN)
    Fh = get(res, "Fh", NaN)

    
    ## Calculate observables
    Fexact = ising_free_energy.(1.0 ./ ts, J)
        
    # relative error in free energy
    re = abs.(F - Fexact) ./ abs.(Fexact) 

    # specific heat and magnetization
    M = -(Fh - F) ./ (h * ts)
    C = SpecificHeat(F,ts)
    
    tp = 0.1:0.01:10.0
    Fexactp = ising_free_energy.(1.0 ./ tp, J)
    Cexact = SpecificHeat(Fexactp,tp)
    Mexact = ising_magnetization.(1.0 ./ tp, J)


    ## Plots

    # create trg images subfolder
    images_path = pwd() * "/images/trg/"
    if !isdir(images_path)
        mkdir(images_path)
    end

    Make2Plot(ts, F, ts, Fexact, Tc, "TRG", "Exact", "Free Energy per site", "T", "F", images_path, "FreeEnergy.png")
    Make2Plot(ts, re, NaN, NaN, Tc, "TRG", "", "Free energy relative error", "T", "ϵ", images_path, "RelativeError.png")
    Make2Plot(ts, C, tp, Cexact, Tc, "TRG", "Exact", "Specific heat per site", "T", "C", images_path, "SpecificHeat.png")
    Make2Plot(ts, M, tp, Mexact, Tc, "TRG", "Exact", "Magnetization", "T", "M", images_path, "Magnetization.png")
    
end