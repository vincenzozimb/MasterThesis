## ---------- Run TRG algorithm ---------- ##

using JLD

include("trg.jl")
include("func.jl")

let 

    ## Physical parameters
    
    # bond dimension (max 15 on my laptop)
    Dcut = 12
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
    
    
    # free energies
    F = - ts .* logZ / 2^(Niter+1);
    Fh = -ts .* logZh / 2^(Niter+1);
    
    
    ## Save results
    
    # add image folder to path
    data_path = pwd() * "/data"
    if !isdir(data_path)
        mkdir(data_path)
    end
    
    # save free energies
    save("data/results.jld", "J", J, "h", h, "ts", ts, "F", F, "Fh", Fh)
    
end