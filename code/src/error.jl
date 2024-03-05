using JLD


include("trg.jl")
include("func.jl")


let

    ## load result from TRG algorithm
    par = load("data/param.jld")

    J = get(par, "J", NaN)
    h = get(par, "h", NaN)

    Dcut = 2:1:15
    Niter = 15 # 11 for 64 x 64 lattice (4096 spins), 13 for  128 x 128 (16384), 15 for 256 x 256 (65536)
    Nspins = 2^(Niter+1)


    ## TRG run
    A = tensor_chess(1/Tc, J, 0.0)
    logZ = []
    for i in Dcut
        y = trg(A, i, Niter)
        push!(logZ, y)
    end 

    F = -(Tc / Nspins) * logZ
    Fexact = ising_free_energy.(1.0 / Tc, J)
    re = abs.(F .- Fexact) ./ abs(Fexact)

    MakePlotLog(Dcut, re, NaN, "", "Relative error at critical point", L"\chi", L"ε", missing, missing)

end