using Plots
using LaTeXStrings

include("trg.jl")

let 

    ## parameters
    Dcut = 4 #max 15
    Niter = 250

    J = 1.0
    Tc = 2.0 / (log(1.0+sqrt(2.0))) # Tc ≈ 2.2691853
    T1 = 0.5 * Tc
    T2 = 2.0 * Tc

    A = tensor_chess(J/T1,0.0)
    tnorm1 = trg_flow(A,Dcut,Niter)
    # A = tensor_chess(J/T2,0.0)
    # above = trg_flow(A,Dcut,Niter)
    
    
    ## plot
    plt = scatter(tnorm1, ms=2, color=:red)
    # scatter!(tnorm1, ms=2, color=:blue)
    
end