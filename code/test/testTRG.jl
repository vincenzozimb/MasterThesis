using LinearAlgebra:svd, Diagonal
using TensorOperations
using Plots
using LaTeXStrings
using QuadGK

let
    
    function Fderivative(F, T)
        delta_T = diff(T)
        dF_dT = diff(F) ./ delta_T
        d2F_dT2 = diff(dF_dT) ./ delta_T[1:end-1]    
        return -T[1:end-2] .* d2F_dT2
    end

    function ising_free_energy(β::Real, J::Real=1.0)
        k = β * J
        c = cosh(2 * k)
        s = sinh(2 * k)
        xmin = 0.0
        xmax = π
        integrand(x) = log(c^2 + √(s^4 + 1 - 2 * s^2 * cos(2.0*x)))
        integral, err = quadgk(integrand, xmin, xmax)::Tuple{Float64,Float64}
        return -(log(2) + integral / π) / (2 * β)
    end

    function TRG(K::Float64, Dcut::Int, no_iter::Int)
        D = 2
        inds = collect(1:D)

        T = zeros(Float64, D, D, D, D)
        M = [[sqrt(cosh(K)) sqrt(sinh(K))];
            [sqrt(cosh(K)) -sqrt(sinh(K))];
            ]
        for i in inds, j in inds, k in inds, l in inds
            for a in inds
                T[i, j, k, l] += M[a, i] * M[a, j] * M[a, k] * M[a, l]
            end
        end

        lnZ = 0.0
        for n in collect(1:no_iter)

            #println(n, " ", maximum(T), " ", minimum(T))
            maxval = maximum(T)
            T = T/maxval
            lnZ += 2^(no_iter-n+1)*log(maxval)

            D_new = min(D^2, Dcut)

            Ma = reshape(permutedims(T, (3, 2, 1, 4)),  (D^2, D^2))
            Mb = reshape(permutedims(T, (4, 3, 2, 1)),  (D^2, D^2))

            F = svd(Ma)

            S1 = reshape(F.U[:,1:D_new]*Diagonal(sqrt.(F.S[1:D_new])), (D, D, D_new))
            S3 = reshape(Diagonal(sqrt.(F.S[1:D_new]))*F.Vt[1:D_new, :], (D_new, D, D))
            F = svd(Mb)
            S2 = reshape(F.U[:,1:D_new]*Diagonal(sqrt.(F.S[1:D_new])), (D, D, D_new))
            S4 = reshape(Diagonal(sqrt.(F.S[1:D_new]))*F.Vt[1:D_new, :], (D_new, D, D))

            @tensor T_new[r, u, l, d] := S1[w, a, r] * S2[a, b, u] * S3[l, b, g] * S4[d, g, w]

            D = D_new
            T = T_new
        end
        trace = 0.0
        for i in collect(1:D)
            trace += T[i, i, i, i]
        end
        lnZ += log(trace)
    end

    Dcut = 10
    n = 12

    ts = 0.1:0.1:5;
    β = inv.(ts);
    J = 1.0;
    @show "=====TRG======"
    logZ = []
    cnt = size(β,1)
    for K in β
        t = J/K
        #T = Ising( K )
        y = TRG(K, Dcut, n);
        #@show lnZ
        # println(1/K, " ", y/2^n)
        push!(logZ,y/2^n)
        print("\r count=$cnt")
        flush(stdout)
        cnt -= 1
    end
    F = - ts.* logZ
    Tc = 2.0 / (log(1.0+sqrt(2.0))) # Tc ≈ 2.2691853


    # pl1 = scatter(ts, F, ms=2, label="TRG")
    # vline!([Tc], line=:red, label=L"T_c")   
    Fexact = ising_free_energy.(1.0 ./ ts, J)
    # plot!(ts, Fexact, label="Exact")
    # title!("Free energy per site")
    # xlabel!("T")
    # ylabel!("F")
    # display(pl1)


    ## relative error in free energy
    re = abs.(F - Fexact) ./ abs.(Fexact)
    
    pl2 = scatter(ts, re, ms=2, label="TRG")
    vline!([Tc], line=:red, label=L"T_c")
    title!("Free energy relative error D=$Dcut")
    xlabel!("T")
    ylabel!("ϵ")
    display(pl2)

    # ## specific heat
    # C = Fderivative(F,ts)
    # tp = 0.1:0.01:5.0
    # Fexact = ising_free_energy.(1.0 ./ tp, J)
    # Cexact = Fderivative(Fexact,tp)
    
    # pl3 = scatter(ts[1:end-2], C, ms=2, label="TRG")
    # vline!([Tc], line=:red, label=L"T_c")
    # plot!(tp[1:end-2], Cexact, lw=:2, label="Exact")
    # title!("Specific heat")
    # xlabel!("T")
    # ylabel!("C")
    # display(pl3)

end