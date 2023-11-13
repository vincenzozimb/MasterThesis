using ITensors
using QuadGK
using ForwardDiff

## Square lattice Ising exact free energy 
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


## Square lattice Ising exact magnetization
function ising_magnetization(β::Real)
    β > βc && return (1 - sinh(2 * β)^(-4))^(1 / 8)
    return 0.0
end

# dlnZ(K::Vector) = ForwardDiff.gradient(K->lnZ(K[1]), K)[1];
# dlnZ2(K::Vector) = ForwardDiff.gradient(K->dlnZ(K), K)[1];


## function that calculates the elementary tensor in the chessboard representation
function tensor_chess(J::Real, h::Real)
    
    # The coupling constants are actually the adimensional ones. That is J is actually βJ and h is actually βh

    dim0 = 2
    s = Index(dim0,"scale=0")
    
    # indeces
    l = addtags(s,"left");
    d = addtags(s,"down");
    r = addtags(s,"right");
    u = addtags(s,"up");

    # lambda function that converts array integers into Ising spin integers (1,2) -> (1,-1)
    Sig = (x) -> 1-2*(x-1)
        
    # define and initialize tensor
    A = ITensor(l,d,r,u)
    
    for sl in 1:dim0
        for sd in 1:dim0
            for sr in 1:dim0
                for su in 1:dim0
                    Eint = Sig(sl)*Sig(sd) + Sig(sd)*Sig(sr) + Sig(sr)*Sig(su) + Sig(su)*Sig(sl)
                    Emag = (Sig(sl) + Sig(sd) + Sig(sr) + Sig(su)) / 2.0
                    P = exp(J*Eint + h*Emag)
                    A[l => sl, d => sd, r => sr, u => su] = P
                end
            end
        end
    end

    return A

end


## function that calculates the log-partition function using TRG
function trg(A::ITensor, maxdim::Int, topscale::Int)
        
    # check A tensor. It is configured to Ising models. It can be modified
    @assert size(A) == (2,2,2,2)
    @assert length(inds(A)) == 4
    @assert hastags(tags(inds(A)[1]), "left,scale=0")
    @assert hastags(tags(inds(A)[2]), "down,scale=0")
    @assert hastags(tags(inds(A)[3]), "right,scale=0")
    @assert hastags(tags(inds(A)[4]), "up,scale=0")

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
    end

    return log(Z)
    
end



