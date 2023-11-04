using ITensors
using QuadGK
using Plots

# function that calculates the partition function using TRG
function trg(T::Float64, maxdim::Int, topscale::Int)
    
    dim0 = 2
    s = Index(dim0,"scale=0")
    
    # indeces
    l = addtags(s,"left");
    r = addtags(s,"right");
    u = addtags(s,"up");
    d = addtags(s,"down");

    # lambda function that converts array integers into Ising spin integers 
    Sig = (x) -> 1-2*(x-1)
        
    # define and initialize tensor
    A = ITensor(l,r,u,d)
    
    for sl in 1:dim0
        for sd in 1:dim0
            for sr in 1:dim0
                for su in 1:dim0
                    E = Sig(sl)*Sig(sd) + Sig(sd)*Sig(sr) + Sig(sr)*Sig(su) + Sig(su)*Sig(sl)
                    P = exp(-E / T)
                    A[l => sl, d => sd, r => sr, u => su] = P
                end
            end
        end
    end

    # TRG algorithm loop
    Z = 1.0
    for scale in 1:topscale
        println("\n---------- Scale $scale -> $(1 + scale)  ----------")
        Fl, Fr = factorize(A, (r,d); maxdim=maxdim, tags="left,scale=$scale")
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
        r = r_new
        u = u_new
        d = d_new

        trace = A*delta(l,r)*delta(u,d)
        A /= trace
        Z *= trace
        # @info Z
    end

    return scalar(Z)

end


# test 

T = 2.0
maxdim = 12
topscale = 10

# Z = trg(T, maxdim, topscale)
# println(Z)


function ising_free_energy(β::Real, J::Real=1.0)
    k = β * J
    c = cosh(2 * k)
    s = sinh(2 * k)
    xmin = 0.0
    xmax = π
    integrand(x) = log(c^2 + √(s^4 + 1 - 2 * s^2 * cos(x)))
    integral, err = quadgk(integrand, xmin, xmax)::Tuple{Float64,Float64}
    return -(log(2) + integral / π) / (2 * β)
end
  
function ising_magnetization(β::Real)
    β > βc && return (1 - sinh(2 * β)^(-4))^(1 / 8)
    return 0.0
end


ts = 0.1:0.01:3.0
@show "=====TRG======"
Z = []
for K in ts
    y = trg(K, maxdim, topscale)
    push!(Z,y)
end
F = - ts .* log.(Z)
beta = 1 ./ ts

pl1 = plot(ts, F)

Fexact = ising_free_energy.(1.0 ./ ts)
pl2 = plot!(ts, Fexact)
display(pl1)

