using ITensors
using QuadGK
using Plots
using LaTeXStrings

let

    ## function that calculates the partition function using TRG
    function trg(K::Float64, maxdim::Int, topscale::Int)
        
        dim0 = 2
        s = Index(dim0,"scale=0")
        
        # indeces
        l = addtags(s,"left");
        d = addtags(s,"down");
        r = addtags(s,"right");
        u = addtags(s,"up");

        # lambda function that converts array integers into Ising spin integers 
        Sig = (x) -> 1-2*(x-1)
            
        # define and initialize tensor
        A = ITensor(l,d,r,u)
        
        for sl in 1:dim0
            for sd in 1:dim0
                for sr in 1:dim0
                    for su in 1:dim0
                        E = Sig(sl)*Sig(sd) + Sig(sd)*Sig(sr) + Sig(sr)*Sig(su) + Sig(su)*Sig(sl)
                        P = exp(K*E)
                        A[l => sl, d => sd, r => sr, u => su] = P
                    end
                end
            end
        end

        # TRG algorithm loop
        Z = 1.0
        for scale in 1:topscale + 1
            # println("\n---------- Scale $scale -> $(1 + scale)  ----------")
            Fl, Fr = factorize(A, (d,r); maxdim=maxdim, tags="left,scale=$scale")
            Fu, Fd = factorize(A, (l,d); maxdim=maxdim, tags="up,scale=$scale")
            
            l_new = commoninds(Fl,Fr)
            r_new = replacetags(l_new, "left", "right")
            Fr *= delta(l_new, r_new)

            u_new = commoninds(Fu,Fd)
            d_new = replacetags(u_new, "up", "down")
            Fd *= delta(u_new, d_new)

            Fl *= delta(r,l)
            Fu *= delta(d,u)
            Fr *= delta(l,r)
            Fd *= delta(u,d)
            
            A = Fl * Fu * Fr * Fd

            l = l_new
            d = d_new
            r = r_new
            u = u_new
            
            trace = A * delta(l,r) * delta(u,d)
            trace = scalar(trace)
            A /= trace
            Z *= trace ^ (1.0 / (2.0^(1+scale)))
            # @info Z
        end

        return Z
        
    end

    
    ## Ising exact free energy 
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
    
    ## heat capacity
    function heat_capacity(F, T)
        delta_T = diff(T)
        dF_dT = diff(F) ./ delta_T
        d2F_dT2 = diff(dF_dT) ./ delta_T[1:end-1]    
        return -T[1:end-2] .* d2F_dT2
    end
    

    ## run
    maxdim = 6
    topscale = 12
    
    Tc = 2.0 / (log(1.0+sqrt(2.0))) # Tc ≈ 2.2691853

    J = 1.0
    ts = 0.1:0.2:10.0
    ks = J ./ ts
    
    @show "=====TRG======"
    pf = []
    cnt = size(ks,1)
    for k in ks
        y = trg(k, maxdim, topscale)
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

    # specific heat
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

end