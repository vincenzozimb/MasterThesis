## ---------- Tensor Renormalization Group algorithm ---------- ## 

using ITensors

include("func.jl")


## elementary tensor in the chessboard representation
function tensor_chess(beta::Real, J::Real, h::Real)
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
                    P = exp(beta*J*Eint + beta*h*Emag)
                    A[l => sl, d => sd, r => sr, u => su] = P
                end
            end
        end
    end
    return A
end


## Elementary tensor in the decomposition W = M M^t
function tensor_eig(K::Real)
    D = 2
    inds = collect(1:D)
    T = zeros(Float64, D, D, D, D)
    M = [[sqrt(cosh(K)) sqrt(sinh(K))];
        [sqrt(cosh(K)) -sqrt(sinh(K))];
        ]
    M = M'
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


## Check tensor properties
function check_tensor(A::ITensor)
    # check A tensor. It is configured to Ising models. It can be modified
    @assert size(A) == (2,2,2,2)
    @assert length(inds(A)) == 4
    @assert hastags(tags(inds(A)[1]), "left,scale=0")
    @assert hastags(tags(inds(A)[2]), "down,scale=0")
    @assert hastags(tags(inds(A)[3]), "right,scale=0")
    @assert hastags(tags(inds(A)[4]), "up,scale=0")
    return true
end


## Calculate the log-partition function using TRG
function trg(A::ITensor, Dcut::Int, Niter::Int)

    # check A tensor
    @assert check_tensor(A)

    # retrieve indeces
    l, d, r, u = filterinds(A)
    @assert hassameinds((l,d,r,u), A)

    # TRG algorithm loop
    lnZ = 0.0
    #initial normalization
    trace = A * delta(l,r) * delta(u,d)
    trace = scalar(trace)
    A /= trace
    lnZ += log(trace) * 2^(Niter)
    for n in 1:Niter
        # tensor decomposition
        Fl, Fr = factorize(A, (d,r); maxdim=Dcut, tags="left,scale=$n")
        Fu, Fd = factorize(A, (l,d); maxdim=Dcut, tags="up,scale=$n")
        
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
        
        # calculate coarse grained tensor
        A = Fl * Fu * Fr * Fd
        
        l = l_new
        d = d_new
        r = r_new
        u = u_new

        # normalize tensor and update
        trace = A * delta(l,r) * delta(u,d)
        trace = scalar(trace)
        A /= trace
        lnZ += log(trace) * 2^(Niter-n)
    end
    return lnZ
end