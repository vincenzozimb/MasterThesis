# using LinearAlgebra
# using TensorOperations
using Plots

include("../src/IsingTRG/trg.jl") 

let 
    
    J = 1.0
    dim0 = 2
    i = Index(dim0,"scale=0")
    l = addtags(i,"left");
    d = addtags(i,"down");
    r = addtags(i,"right");
    u = addtags(i,"up");
        
    # define and initialize tensor
    M = [sqrt(cosh(J)) sqrt(sinh(J)); sqrt(cosh(J)) -sqrt(sinh(J))]
    M1 = ITensor(M,i,l)
    M2 = ITensor(M,i,d)
    M3 = ITensor(M,i,r)
    M4 = ITensor(M,i,u)

    @info M * M'
    @info inds(M1 * M2 * M3 * M4)



    # T = 1.0
    # Tc = 2.0 / (log(1.0+sqrt(2.0))) # Tc ≈ 2.2691853
    # Tl = Tc - 0.1
    # Th = Tc + 0.1

    # A = tensor_chess(1/T, 0.0)
    # Al = tensor_chess(1/Tl, 0.0)
    # Ac = tensor_chess(1/Tc, 0.0)
    # Ah = tensor_chess(1/Th, 0.0)

    # niter = 150
    # lnZ, n = trg(A, 6, niter)
    # lnZl, nl = trg(Al, 6, niter)
    # lnZc, nc = trg(Ac, 6, niter)
    # lnZh, nh = trg(Ah, 6, niter)
    
    # plot(n, label="T", marker=:square)
    # plot!(nl, label="T low", marker=:circle)
    # plot!(nc, label="T critical", marker=:square)
    # plot!(nh, label="T high")

end


    
