## ---------- Run TRG algorithm ---------- ##

using JLD

include("trg.jl")
include("func.jl")

let 

    ## Physical parameters
    
    # bond dimension (max 15 on my laptop)
    Dcut = 15
    # number of iterations
    Niter = 15 # 11 for 64 x 64 lattice (4096 spins), 13 for  128 x 128 (16384), 15 for 256 x 256 (65536)
    # coupling constant
    J = 1.0
    # external field
    h = 0.01


    ## Simulation parameters

    # temperature range
    ts = 0.5:0.15:10.0
    # βJ range
    bs = 1.0 ./ ts
    

    ## TRG run
    welcome = "=====TRG======" 
    @show welcome
    logZ = []
    logZh = []
    cnt = size(bs,1)
    for b in bs
        A = tensor_chess(b, J, 0.0)
        y = trg(A, Dcut, Niter)
        push!(logZ,y)
        A = tensor_chess(b, J, h)
        y = trg(A, Dcut, Niter)
        push!(logZh,y)
        print("\r count=$cnt")
        flush(stdout)
        cnt -= 1
    end
    println()
    
    
    # free energies per spin
    F = - ts .* logZ / 2^(Niter+1);
    Fh = -ts .* logZh / 2^(Niter+1);
    
    
    ## Save results
    
    # add image folder to path
    data_path = pwd() * "/data"
    if !isdir(data_path)
        mkdir(data_path)
    end
    
    # save parameters and results
    save("data/param.jld", "Niter", Niter, "J", J, "h", h)
    save("data/resultsTrg.jld", "ts", ts, "F", F, "Fh", Fh, "logZ", logZ, "logZh", logZh)
    
end