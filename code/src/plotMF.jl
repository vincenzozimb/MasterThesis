using JLD
using Statistics

include("func.jl")

let 

    # ## Load data
    # par = load("data/param.jld")
    # mf = load("data/MeanField.jld")

    # idx = get(mf, "idx", NaN)
    # ts = get(mf, "ts", NaN)
    # magn_sample = get(mf, "magn_sample", NaN)
    # ener_sample = get(mf, "ener_sample", NaN)
    # nsamples = get(mf, "nsamples", NaN)

    # Niter = get(par, "Niter", NaN)
    # Nspins = 2^(Niter+1)
    # @assert isinteger(sqrt(Nspins))


    # ## Calculate mean and standart deviation (samples in the form sample[nsamples,ts]
    # M = mean(magn_sample) ./ Nspins
    # dM = std(magn_sample) ./ Nspins

    # E = mean(ener_sample) ./ Nspins
    # VarE = var(ener_sample) ./ Nspins^2
    # dE = sqrt.(VarE)
   
    
    # ## Specific heat as numerical derivative
    # C = diff(E) ./ diff(ts[idx])
    # append!(C,NaN) 
    
    # # error calculation
    # a1 = vec(mean(ener_sample[:,2:end] .* ener_sample[:,1:end-1], dims=1))
    # b1 = vec(mean(ener_sample[:,2:end], dims=1)) .* vec(mean(ener_sample[:,1:end-1], dims=1))
    
    # a2 = vec(mean(ener_sample[:,3:end] .* ener_sample[:,1:end-2], dims=1))
    # b2 = vec(mean(ener_sample[:,3:end], dims=1)) .* vec(mean(ener_sample[:,1:end-2], dims=1))
    
    # G1 = a1 .- b1
    # append!(G1, NaN)
    # G2 = a2 .- b2
    # append!(G2, [NaN, NaN])
     
    # dC = -(G2 .- 2*G1 .+ VarE) ./ (diff(ts[idx])[1] ^2)
    # dC = sqrt.(abs.(dC)) ./ Nspins
    
    # # # alternative with the error propagation
    # # dC = sqrt.(VarE[2:end] .+ VarE[1:end-1]) ./ diff(ts[idx])


    # ## Specific heat as energy fluctuations (corect formula but no correct result), and bootstrap needed for error
    # C2 = Nspins * (VarE ./ (ts[idx] .^2))
    
    
    # ## Plots

    # # create MF images subfolder
    # images_path = pwd() * "/images/MF"
    # if !isdir(images_path)
    #     mkdir(images_path)
    # end

    # MakeErrorPlot(ts[idx], M, dM, 4, "Mean field", "Magnetization with mean field", "T", "M", "MF/magnetization.png")
    # MakeErrorPlot(ts[idx], E, dE, 4, "Mean field", "Energy per site with mean field", "T", "E", "MF/energy.png")
    # MakeErrorPlot(ts[idx], C, dC, 4, "Mean field", "Specific heat per site with num der", "T", "E", "MF/specificheat.png")
    # MakePlot(ts[idx], C2, NaN, NaN, 4, "Mean field", "", "Specific heat per site with mean field", "T", "C", "MF/specificheatFluct.png")
       
end