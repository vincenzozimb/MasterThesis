# using DSP
# using Plots
# using ITensors

# include("../src/IsingTRG/trg.jl") 

using LinearAlgebra
using TensorOperations

function TRG(K::Float64,iteration::Int64,maxD::Int64)

    T = zeros(Float64,2,2,2,2)

    no2s(no::Int64) = -2*no+3
    sumS(s0::Int64,s1::Int64,s2::Int64,s3::Int64) = no2s(s0)*no2s(s1)+no2s(s1)*no2s(s2)+no2s(s2)*no2s(s3)+no2s(s3)*no2s(s0)

    for s0 = 1:2,s1 = 1:2,s2 = 1:2,s3 = 1:2
        @inbounds T[s3,s2,s1,s0] = exp(K*sumS(s3,s2,s1,s0))
    end

    D = 2
    epsilon = 0.0

    checkBigger(M::Float64) = M>=epsilon

    lnZ = 0.0

    @inbounds for n = 1:iteration
        maxVal = maximum(T)
        T /= maxVal
        lnZ += 2^(iteration-n+1)*log(maxVal)
        TypeA = reshape(T,(D^2,D^2))
        TypeB = reshape(permutedims(T,[1 4 2 3]),(D^2,D^2))
        ua,sa,va = svd(TypeA)
        ub,sb,vb = svd(TypeB)
        cut = min(count(checkBigger,sb),count(checkBigger,sa),maxD)
        S2 = reshape(view(ub,:,1:cut)*Diagonal(sqrt.(view(sb,1:cut))),(D,D,cut))
        S1 = reshape(Diagonal(sqrt.(view(sb,1:cut)))*view(transpose(vb),1:cut,:),(cut,D,D))
        S3 = reshape(view(ua,:,1:cut)*Diagonal(sqrt.(view(sa,1:cut))),(D,D,cut))
        S4 = reshape(Diagonal(sqrt.(view(sa,1:cut)))*view(transpose(va),1:cut,:),(cut,D,D))

        @tensor T[a,e,g,n] := S1[a,b,c]*S3[c,d,e]*S2[f,d,g]*S4[n,f,b]
        #T = Tp
        D = cut
    end

    trace = 0.0
    for i = 1:D
        @inbounds trace += T[i,i,i,i]
    end
    lnZ += log(trace)

    println("K:",K)
    println("lnZ/N:",lnZ/(2^(iteration+1)))
end

function main(args)
    if isempty(args)
        K = 1.0
        iters = 20
        cut = 20
    else
        K = parse(Float64,args[1])
        iters = parse(Int64,args[2])
        cut = parse(Int64,args[3])
    end
    TRG(K,iters,cut)
end

main(ARGS)



# let 
    
#     beta = 1.0
#     J = 1.0 * beta
#     h = 0.1 * beta

#     A = tensor_chess(J, h)
    
#     @assert size(A) == (2,2,2,2)
#     @assert length(inds(A)) == 4
#     @assert hastags(tags(inds(A)[1]), "left,scale=0")
#     @assert hastags(tags(inds(A)[2]), "down,scale=0")
#     @assert hastags(tags(inds(A)[3]), "right,scale=0")
#     @assert hastags(tags(inds(A)[4]), "up,scale=0")


    
    

    
# end




    # # Generate a noisy signal
    # sample_rate = 1000.0  # Sample rate in Hz
    # t = 0:1/sample_rate:1  # Time vector
    # x = sin.(2π * 100 * t) + 0.5 * randn(length(t))  # Noisy signal (sinusoid + noise)

    # # Plot the original signal
    # pl = plot(t, x, label="Noisy Signal")
    # # plot!(t, filtered_signal, label="Filtered Signal")
    # xlabel!("Time (s)")
    # ylabel!("Amplitude")
    # # title!("Low-Pass Filtered Signal")
    
    # # filter
    # f_cutoff = 0.2
    # responsetype = Lowpass(10; fs=sample_rate)
    # designmethod = Butterworth(4)
    # xf = filt(digitalfilter(responsetype, designmethod), x)  

    # # Plot the filtered signal
    # plot!(t, xf, label="Filtered Signal")
    # title!("Low-Pass Filtered Signal")
    # display(pl)

    
