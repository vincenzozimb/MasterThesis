using JLD
using Statistics

include("func.jl")

let 

    ## load result from TRG algorithm
    par = load("data/param.jld")
    mf = load("data/MeanField.jld")

    idx = get(mf, "idx", NaN)
    ts = get(mf, "ts", NaN)
    magn_sample = get(mf, "magn_sample", NaN)
    energy_sample = get(mf, "energy_sample", NaN)
    nsamples = get(mf, "nsamples", NaN)

    Niter = get(par, "Niter", NaN)
    Nspins = 2^(Niter+1)
    @assert isinteger(sqrt(Nspins))


    ## calculate mean and standart deviation
    magn = vec(mean(magn_sample, dims=1)) ./ Nspins
    dmagn = vec(std(magn_sample, dims=1)) ./ Nspins # to check

    energy = vec(mean(energy_sample, dims=1))
    energyVar = vec(var(energy_sample, dims=1)) 
    specificheat = (energyVar ./ (ts[idx] .^2)) ./ Nspins
    energy ./= Nspins


    ## Plots

    # create trg images subfolder
    images_path = pwd() * "/images/MF"
    if !isdir(images_path)
        mkdir(images_path)
    end

    x = 0:0.1:9
    y = fill(1.5,length(0:0.1:9))

    MakePlot(ts[idx], magn, NaN, NaN, 4, "Mean field", "", "Magnetization with mean field", "T", "M", "MF/magnetization.png")
    MakePlot(ts[idx], energy, NaN, NaN, 4, "Mean field", "", "Energy per site with mean field", "T", "E", "MF/energy.png")
    MakePlot(ts[idx], specificheat, x, y, 4, "Mean field", "", "Specific heat per site with mean field", "T", "C", "MF/specificheat.png")
    
end