using Plots

include("../src/IsingTRG/trg.jl")

let 

    ## parameters
    Dcut = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
    Dplot = map(string, Dcut)
    Niter = 15

    J = 1.0
    Tc = 2.0 / (log(1.0+sqrt(2.0))) # Tc ≈ 2.2691853
    T1 = 0.5 * Tc
    T2 = 2.0 * Tc

    
    B = tensor_chess(J/T1,0.0)
    t = Array{Float64}(undef, 0)
    
    for i in eachindex(Dcut)
        value = @elapsed trg(B, Dcut[i], Niter)
        push!(t, value)
    end
    
    plt = scatter(Dplot, t, ms=6, color=:red, label=nothing)
    # plot!(Dcut, Dcut.^6)
    xlabel!("Dcut")
    ylabel!("time for trg run [s]")
    display(plt)

    images_path = pwd() * "/images"
    if !isdir(images_path)
        mkdir(images_path)
    end
    images_path = pwd() * "/images"
    savefig(joinpath(images_path, "TIMEvsD.png"))
    
    
end
