using Plots
using LaTeXStrings

include("trg.jl")

let 

    ## parameters
    Dcut = 15 #max 15
    Niter = 15

    J = 1.0
    Tc = 2.0 / (log(1.0+sqrt(2.0))) # Tc ≈ 2.2691853
    T1 = 0.5 * Tc
    T2 = 2.0 * Tc

    A = tensor_chess(J/T1,0.0)
    @time trg(A, Dcut, Niter)
    nothing
    
    


end