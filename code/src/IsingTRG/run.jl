using Plots
using LaTeXStrings

include("trg.jl")

let   

    ## functions
    function Fderivative(F, T)
        delta_T = diff(T)
        dF_dT = diff(F) ./ delta_T
        d2F_dT2 = diff(dF_dT) ./ delta_T[1:end-1]    
        return -T[1:end-2] .* d2F_dT2
    end

    function moving_average(data, window_size)
        n = length(data)
        result = zeros(n)
    
        for i in 1:n
            lower = max(1, i - window_size ÷ 2)
            upper = min(n, i + window_size ÷ 2)
            result[i] = sum(data[lower:upper]) / (upper - lower + 1)
        end
    
        return result
    end

    ## parameters
    maxdim = 6 #max 15
    topscale = 12
    
    Tc = 2.0 / (log(1.0+sqrt(2.0))) # Tc ≈ 2.2691853

    J = 1.0
    h = 0.01
    ts = 0.1:0.15:10.0
    ks = J ./ ts
    hs = 0.0:0.02:0.1

    
    @show "=====TRG======"
    logZ = []
    logZh = []
    cnt = size(ks,1)
    for k in ks
        A = tensor_chess(k, 0.0)
        y = trg(A, maxdim, topscale)
        push!(logZ,y)
        A = tensor_chess(k, h)
        y = trg(A, maxdim, topscale)
        push!(logZh,y)
        print("\r count=$cnt")
        flush(stdout)
        cnt -= 1
    end
    println()
    F = - ts .* logZ
    Fh = -ts .* logZh

    # #filter by moving average
    # window_size = 5
    # F = moving_average(F, window_size)
    # Fh = moving_average(Fh, window_size)
    # F[1:2] .= NaN
    # Fh[1:2] .= NaN
    # F[end-1:end] .= NaN
    # Fh[end-1:end] .= NaN
    
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
    # re = (Fexact - F) ./ Fexact
    
    pl2 = scatter(ts, re, ms=2, label="TRG")
    vline!([Tc], line=:red, label=L"T_c")
    title!("Free energy relative error")
    xlabel!("T")
    ylabel!("ϵ")
    display(pl2)

    ## specific heat
    C = Fderivative(F,ts)
    tp = 0.1:0.01:10.0
    Fexact = ising_free_energy.(1.0 ./ tp, J)
    Cexact = Fderivative(Fexact,tp)

    pl3 = scatter(ts[1:end-2], C, ms=2, label="TRG")
    vline!([Tc], line=:red, label=L"T_c")
    plot!(tp[1:end-2], Cexact, lw=:2, label="Exact")
    title!("Specific heat")
    xlabel!("T")
    ylabel!("C")
    display(pl3)


    ## magnetization
    M = -(Fh - F) ./ (h * ts)
    tp = 0.1:0.01:10.0
    Fexact = ising_free_energy.(1.0 ./ tp, J)
    Mexact = ising_magnetization.(1.0 ./ tp)

    pl4 = scatter(ts, M, ms=2, label="TRG")
    vline!([Tc], line=:red, label=L"T_c")
    plot!(tp, Mexact, lw=:2, label="Exact")
    title!("Magnetization")
    xlabel!("T")
    ylabel!("M")
    display(pl4)


    # ## relative error in the specific heat
    # Fexact = ising_free_energy.(1.0 ./ ts, J)
    # Cexact = Fderivative(Fexact,ts)

    # re = abs.(C - Cexact) ./ abs.(Cexact)
    
    # pl4 = scatter(ts, re, ms=2, label="TRG")
    # vline!([Tc], line=:red, label=L"T_c")
    # title!("Specific heat relative error")
    # xlabel!("T")
    # ylabel!("ϵ")
    # display(pl4)
    
end