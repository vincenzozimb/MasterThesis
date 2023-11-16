using Plots
using ITensors
using LaTeXStrings

let 
    
    ## ising free energy
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

    ## function that calculates the elementary tensor with the W = M * M^t, without magnetic field for now
    function tensor_eig(K::Float64)
        D = 2
        inds = collect(1:D)
        T = zeros(Float64, D, D, D, D)
        M = [[sqrt(cosh(K)) sqrt(sinh(K))];
            [sqrt(cosh(K)) -sqrt(sinh(K))];
            ]
        for i in inds, j in inds, k in inds, m in inds
            for a in inds
                T[i, j, k, m] += M[a, i] * M[a, j] * M[a, k] * M[a, m]
            end
        end
        s = Index(D,"scale=0")
        l = addtags(s,"left");
        d = addtags(s,"down");
        r = addtags(s,"right");
        u = addtags(s,"up");
        return ITensor(T,l,d,r,u)
    end

    # ## function that perform Tensor Renormalization Group to the tensor given
    # function trg(A::ITensor, chi::Int, n_iter::Int)
    #     # check A tensor. It is configured to Ising models. It can be modified
    #     @assert size(A) == (2,2,2,2)
    #     @assert length(inds(A)) == 4
    #     @assert hastags(tags(inds(A)[1]), "left,scale=0")
    #     @assert hastags(tags(inds(A)[2]), "down,scale=0")
    #     @assert hastags(tags(inds(A)[3]), "right,scale=0")
    #     @assert hastags(tags(inds(A)[4]), "up,scale=0")
    #     # retrieve indeces
    #     l, d, r, u = inds(A)
    #     @assert hassameinds((l,d,r,u), A)
    #     # TRG algorithm loop
    #     lnZ = 0.0
    #     for n in 1:n_iter
    #         # renormalize tensor
    #         trace = A * delta(r,l) * delta(u,d)
    #         A = A / trace
    #         lnZ += log(trace) * 2^(n_iter+n-1)
    #         # factorize tensor
    #         Fl, Fr = factorize(A, (d,r); maxdim=maxdim, tags="left,scale=$scale")
    #         Fu, Fd = factorize(A, (l,d); maxdim=maxdim, tags="up,scale=$scale")
    #         # new indeces            
    #         l_new = commoninds(Fl,Fr)
    #         r_new = replacetags(l_new, "left", "right")
    #         Fr *= delta(l_new, r_new)
    #         u_new = commoninds(Fu,Fd)
    #         d_new = replacetags(u_new, "up", "down")
    #         Fd *= delta(u_new, d_new)
    #         # use the correct indeces
    #         Fl *= delta(r,l)
    #         Fu *= delta(d,u)
    #         Fr *= delta(l,r)
    #         Fd *= delta(u,d)
    #         # build coarse-grained tensor            
    #         A = Fl * Fu * Fr * Fd
    #         # update indeces
    #         l = l_new
    #         d = d_new
    #         r = r_new
    #         u = u_new
    #     end
        

    # end

    function trg(T::ITensor; χmax::Int, nsteps::Int, cutoff=0.0, svd_alg="divide_and_conquer")
        sₕ, sᵥ = filterinds(T; plev=0)
        @assert hassameinds((sₕ, sₕ', sᵥ, sᵥ'), T)
      
        # Keep track of the partition function per site
        κ = 1.0
        for n in 1:nsteps
          Fₕ, Fₕ′ = factorize(
            T, (sₕ', sᵥ'); ortho="none", maxdim=χmax, cutoff, tags=tags(sₕ), svd_alg
          )
      
          s̃ₕ = commonind(Fₕ, Fₕ′)
          Fₕ′ *= δ(dag(s̃ₕ), s̃ₕ')
      
          Fᵥ, Fᵥ′ = factorize(
            T, (sₕ, sᵥ'); ortho="none", maxdim=χmax, cutoff, tags=tags(sᵥ), svd_alg
          )
      
          s̃ᵥ = commonind(Fᵥ, Fᵥ′)
          Fᵥ′ *= δ(dag(s̃ᵥ), s̃ᵥ')
      
          T =
            (Fₕ * δ(dag(sₕ'), sₕ)) *
            (Fᵥ * δ(dag(sᵥ'), sᵥ)) *
            (Fₕ′ * δ(dag(sₕ), sₕ')) *
            (Fᵥ′ * δ(dag(sᵥ), sᵥ'))
      
          sₕ, sᵥ = s̃ₕ, s̃ᵥ
      
          trT = abs((T * δ(sₕ, sₕ') * δ(sᵥ, sᵥ'))[])
          T = T / trT
          κ *= trT^(1 / 2^n)
        end
        return κ
    end

    function tensor_ising(K::Float64)
        dim0 = 2
        sₕ = Index(dim0) 
        sᵥ = Index(dim0)
    
        A = tensor_eig(K)
        A = Array(A,inds(A))
        A = ITensor(A,sₕ,sₕ',sᵥ,sᵥ')
        return A
    end


    Dcut = 4
    n = 30

    # A = tensor_ising(1.0)
    # y = trg(A;χmax=Dcut,nsteps=n)

    ts = 0.5:0.1:5;
    β = inv.(ts);
    J = 1.0;
    @show "=====TRG======"
    logZ = []
    cnt = size(β,1)
    for K in β
        t = J*K
        A = tensor_ising(t)
        y = trg(A; χmax=Dcut, nsteps=n)
        #@show lnZ
        # println(1/K, " ", y/2^n)
        # push!(logZ,y/2^n)
        push!(logZ,y)
        print("\r count=$cnt")
        flush(stdout)
        cnt -= 1
    end
    F = - ts.* logZ
    Tc = 2.0 / (log(1.0+sqrt(2.0))) # Tc ≈ 2.2691853


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
    C = Fderivative(F,ts)
    tp = 0.1:0.01:5.0
    Fexact = ising_free_energy.(1.0 ./ tp, J)
    Cexact = Fderivative(Fexact,tp)
    
    pl3 = scatter(ts[1:end-2], C, ms=2, label="TRG")
    vline!([Tc], line=:red, label=L"T_c")
    plot!(tp[1:end-2], Cexact, lw=:2, label="Exact")
    title!("Specific heat")
    xlabel!("T")
    ylabel!("C")
    display(pl3)









end


    
