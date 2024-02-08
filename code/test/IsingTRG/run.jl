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
    Dcut = 6 #max 15
    Niter = 15
    
    Tc = 2.0 / (log(1.0+sqrt(2.0))) # Tc ≈ 2.2691853

    J = 1.0
    h = 0.01
    ts = 0.1:0.15:10.0
    ks = J ./ ts
    hs = 0.0:0.02:0.1

    
    @show "=====TRG======"
    logZ = []
    logZh = []
    tensors = []
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
    F = - ts .* logZ / 2^(Niter+1)
    Fh = -ts .* logZh / 2^(Niter+1)

    # #filter by moving average
    # window_size = 3
    # F = moving_average(F, window_size)
    # Fh = moving_average(Fh, window_size)
    # F[1] = NaN
    # Fh[1] = NaN
    # F[end] = NaN
    # Fh[end] = NaN
    
    # add images folder to path
    images_path = pwd() * "/images"
    if !isdir(images_path)
        mkdir(images_path)
    end
    
    pl1 = scatter(ts, F, ms=2, label="TRG")
    vline!([Tc], line=:red, label=L"T_c")   
    Fexact = ising_free_energy.(1.0 ./ ts, J)
    plot!(ts, Fexact, label="Exact")
    title!("Free energy per site D=$Dcut Niter=$Niter")
    xlabel!("T")
    ylabel!("F")
    savefig(joinpath(images_path, "FD$Dcut"*"N$Niter.png"))

    ## relative error in free energy
    re = abs.(F - Fexact) ./ abs.(Fexact)
    # re = (Fexact - F) ./ Fexact
    
    pl2 = scatter(ts, re, ms=2, label="TRG")
    vline!([Tc], line=:red, label=L"T_c")
    title!("Free energy relative error D=$Dcut Niter=$Niter")
    xlabel!("T")
    ylabel!("ϵ")
    savefig(joinpath(images_path, "ReD$Dcut"*"N$Niter.png"))

    ## specific heat
    C = Fderivative(F,ts)
    tp = 0.1:0.01:10.0
    Fexact = ising_free_energy.(1.0 ./ tp, J)
    Cexact = Fderivative(Fexact,tp)

    pl3 = scatter(ts[1:end-2], C, ms=2, label="TRG")
    vline!([Tc], line=:red, label=L"T_c")
    plot!(tp[1:end-2], Cexact, lw=:2, label="Exact")
    title!("Specific heat D=$Dcut Niter=$Niter")
    xlabel!("T")
    ylabel!("C")
    savefig(joinpath(images_path, "CD$Dcut"*"N$Niter.png"))


    ## magnetization
    M = -(Fh - F) ./ (h * ts)
    tp = 0.1:0.01:10.0
    Fexact = ising_free_energy.(1.0 ./ tp, J)
    Mexact = ising_magnetization.(1.0 ./ tp)

    pl4 = scatter(ts, M, ms=2, label="TRG")
    vline!([Tc], line=:red, label=L"T_c")
    plot!(tp, Mexact, lw=:2, label="Exact")
    title!("Magnetization D=$Dcut Niter=$Niter")
    xlabel!("T")
    ylabel!("M")
    savefig(joinpath(images_path, "MD$Dcut"*"N$Niter.png"))
       

end