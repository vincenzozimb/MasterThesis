## ---------- MEAN FIELD THEORY FOR SQUARE LATTICE ISING MODEL ---------- ##

using JLD
using Distributions

include("func.jl")


## function for self consistent equation 
function mf(L::Real, β::Real; h::Real=0.0, J::Real=1.0, q::Int=4)
    return L - tanh(β * (h + q*J*L) )
end


## derivative of the function for self consistent equation
function mfder(L::Real, β::Real; h::Real=0.0, J::Real=1.0, q::Int=4)
    return 1.0 - β*J*q / (cosh( β * (h + q*J*L) )^2)
end


## sample configuration from the mean field distribution
function sampleMF(Nspins::Integer, pUp)
    n = Int(sqrt(Nspins))
    spins = hcat([rand(Bernoulli(p), Nspins) for p in pUp]...)
    spins = (2 * spins) .- 1
    spins = reshape(spins, n, n, length(pUp)) 
    return spins
end


## Ising energy function
function IsingEnergy(spins::Array{Int, 3}; J::Float64=1.0, h::Float64=0.0)

    Lx, Ly, nt = size(spins)

    energy_t = zeros(Float64, nt)

    for t in 1:nt
        energy = 0

        for i in 1:Lx
            for j in 1:Ly
                s = spins[i, j, t]
                nb_sum = spins[mod1(i+1, Lx), j, t] +
                         spins[mod1(i-1, Lx), j, t] +
                         spins[i, mod1(j+1, Ly), t] +
                         spins[i, mod1(j-1, Ly), t]

                energy += -J * s * nb_sum - h * s
            end
        end

        energy_t[t] = energy / 2  # divide by 2 to avoid double counting
    end

    return energy_t

end


let 
    
    ## load data
    par = load("data/param.jld")
    
    J = get(par, "J", NaN)
    h = get(par, "h", NaN)
    Niter = get(par, "Niter", NaN)


    ## total number of spins (lattice of sqrt(Nspins) x sqrt(Nspins) spins )
    Nspins = 2^(Niter+1)
    @assert isinteger(sqrt(Nspins))
   

    ## temperature range
    ts = 0.5:0.05:8

    ## calculate zeros (NaN at the critical point)
    x0 = 1.0
    Tcmf = 4.0 # mean field critical temperature
    
    L0 = []
    for t in ts
        t != Tcmf ? y = NewtonRaphson(mf, x0, mfder, (1.0/t,)) : y = NaN 
        push!(L0,y)
    end


    ## calculate mean field distribution (probability is NaN at the critical point)
    q = 4
    heff = (J*q*L0 .+ h) ./ ts
    pUp = exp.(heff) ./ (exp.(heff) .+ exp.(-heff))
    

    ## exclude critical temperature and take a simple sample
    idx = .!isnan.(pUp)
    pUp = pUp[idx]
    spins = sampleMF(Nspins, pUp)
    

    # # calculate and plot sample magnetization and energy
    # images_path = pwd() * "/images/MF"
    # if !isdir(images_path)
    #     mkdir(images_path)
    # end
    # M = vec(sum(spins, dims=(1,2)) ./ Nspins)    
    # MakePlot(ts[idx], M, NaN, NaN, 4, "Mean field", "", "Magnetization with mean field", "T", "M", "MFmagnetizationSample.png")
    # energy = IsingEnergy(spins) ./ Nspins
    # MakePlot(ts[idx], energy, NaN, NaN, 4, "Mean field", "", "Energy per spin with mean field", "T", "E", "MFenergySample.png")


    ## sample from the mean field distribution
    nsamples = 1000
    magn_sample = zeros(Float64,nsamples, length(ts[idx]))
    energy_sample = zeros(Float64,nsamples, length(ts[idx]))
    for i in 1:nsamples
        configuration = sampleMF(Nspins,pUp)
        magn_sample[i, :] = vec(sum(configuration, dims=(1,2)))
        energy_sample[i, :] = IsingEnergy(configuration)
    end


    # add image folder to path
    data_path = pwd() * "/data"
    if !isdir(data_path)
        mkdir(data_path)
    end


    # save mean field sample
    save("data/MeanField.jld", "ts", ts, "idx", idx, "nsamples", nsamples, "magn_sample", magn_sample, "energy_sample", energy_sample)
    

end