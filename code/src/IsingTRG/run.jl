using Plots
using LaTeXStrings

include("trg.jl")
include("util.jl")

let 

    ## functions
    function heat_capacity(F, T)
        delta_T = diff(T)
        dF_dT = diff(F) ./ delta_T
        d2F_dT2 = diff(dF_dT) ./ delta_T[1:end-1]    
        return -T[1:end-2] .* d2F_dT2
    end


    ## parameters
    maxdim = 6 #max 15
    topscale = 12
    
    Tc = 2.0 / (log(1.0+sqrt(2.0))) # Tc ≈ 2.2691853

    J = 1.0
    h = 0.0
    ts = 0.1:0.2:10.0
    ks = J ./ ts
    # hs = (0.0:0.02:0.1) ./ ts
    
    @show "=====TRG======"
    pf = []
    cnt = size(ks,1)
    for k in ks
        y = trg(k, maxdim, topscale, h)
        push!(pf,y)
        print("\r count=$cnt")
        flush(stdout)
        cnt -= 1
    end
    println()
    F = - ts .* log.(pf)
    
    pl1 = scatter(ts, F, ms=2, label="TRG")
    vline!([Tc], line=:red, label=L"T_c")   
    Fexact = ising_free_energy.(1.0 ./ ts, J)
    plot!(ts, Fexact, label="Exact")
    title!("Free energy per site")
    xlabel!("T")
    ylabel!("F")
    display(pl1)
    
    ## relative error in free energy
    re = abs.(F - Fexact) ./ abs.(Fexact)
    
    pl2 = scatter(ts, re, ms=2, label="TRG")
    vline!([Tc], line=:red, label=L"T_c")
    title!("Free energy relative error")
    xlabel!("T")
    ylabel!("ϵ")
    display(pl2)

    ## specific heat
    C = heat_capacity(F,ts)
    tp = 0.1:0.01:10.0
    Fexact = ising_free_energy.(1.0 ./ tp, J)
    Cexact = heat_capacity(Fexact,tp)

    pl3 = scatter(ts[1:end-2], C, ms=2, label="TRG")
    vline!([Tc], line=:red, label=L"T_c")
    plot!(tp[1:end-2], Cexact, lw=:2, label="Exact")
    title!("Specific heat")
    xlabel!("T")
    ylabel!("C")
    display(pl3)

    ## relative error in the specific heat
    Fexact = ising_free_energy.(1.0 ./ ts, J)
    Cexact = heat_capacity(Fexact,ts)

    re = abs.(C - Cexact) ./ abs.(Cexact)
    
    pl4 = scatter(ts, re, ms=2, label="TRG")
    vline!([Tc], line=:red, label=L"T_c")
    title!("Specific heat relative error")
    xlabel!("T")
    ylabel!("ϵ")
    display(pl4)
    
end