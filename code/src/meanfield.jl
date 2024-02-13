## ---------- MEAN FIELD THEORY FOR SQUARE LATTICE ISING MODEL ---------- ##

using JLD
using Distributions
# using Plots

include("func.jl")


# ## function for self consistent equation 
# function mf(L::Real, β::Real, J::Real=1.0, h::Real=0.0)
#     return L - tanh(β * (h + q*J*L) )
# end


# ## derivative of the function for self consistent equation
# function mfder(L::Real, β::Real, J::Real=1.0, h::Real=0.0)
#     return 1.0 - β*J*q / (cosh( β * (h + q*J*L) )^2)
# end


# ## sample configuration from the mean field distribution
# function sampleMF(Nspins::Integer, p::Real)
#     n = Int(sqrt(Nspins))
#     spins = rand(Bernoulli(p), Nspins)
#     spins = (2 * spins) .- 1
#     spins = reshape(spins, n, n) 
#     return spins
# end


# ## Ising energy functions
# function IsingEnergy(spins::Array{Int, 2}; J::Float64=1.0, h::Float64=0.0)
#     Lx, Ly = size(spins)
#     energy = 0
#     for i in 1:Lx
#         for j in 1:Ly
#             s = spins[i, j]
#             nb_sum = spins[mod1(i+1, Lx), j] +
#                         spins[mod1(i-1, Lx), j] +
#                         spins[i, mod1(j+1, Ly)] +
#                         spins[i, mod1(j-1, Ly)]

#             energy += -J * s * nb_sum - h * s
#         end
#     end
#     energy /= 2  # divide by 2 to avoid double counting
#     return energy
# end


# ## Mean Field magnetization
# function MFmagnetization(ts, J=1.0, h=0.0)
#     # (NaN at the critical point)
#     Tcmf = 4.0
#     x0 = 1.0
#     L0 = []
#     for t in ts
#         t != Tcmf ? y = NewtonRaphson(mf, x0, mfder, (1.0/t, J, h)) : y = NaN 
#         push!(L0,y)
#     end
#     return L0
# end



let 
    
    # ## load data
    # par = load("data/param.jld")
    
    # J = get(par, "J", NaN)
    # h = get(par, "h", NaN)
    # Niter = get(par, "Niter", NaN)


    # ## total number of spins (lattice of sqrt(Nspins) x sqrt(Nspins) spins )
    # Nspins = 2^(Niter+1)
    # @assert isinteger(sqrt(Nspins))
   

    # ## temperature range
    # ts = 0.005:0.1:8.0

    # ## calculate magnetization solving the self consistent equation (NaN at the critical point)   
    # L0 = MFmagnetization(ts)


    # ## calculate mean field distribution (probability is NaN at the critical point)
    # heff = (J*q*L0 .+ h) ./ ts
    # pUp = exp.(heff) ./ (exp.(heff) .+ exp.(-heff))
    

    # ## exclude critical temperature and take a simple sample
    # idx = .!isnan.(pUp)
    # pUp = pUp[idx]
    # ts = ts[idx]


    # ## sampling
    # nsamples = 10
    # magn_sample = []
    # ener_sample = []
    # for i in 1:nsamples
    #     spins = sampleMF.(Nspins, pUp)
    #     push!(magn_sample, sum.(spins) ./ Nspins)
    #     push!(ener_sample, IsingEnergy.(spins, J=J, h=h) ./ Nspins)
    # end

    # ener = mean(ener_sample)
    # dener = std(ener_sample)
    # magn = mean(magn_sample)
    # dmagn = std(magn_sample)

    # # add image folder to path
    # data_path = pwd() * "/images/MF"
    # if !isdir(data_path)
    #     mkdir(data_path)
    # end

    # MakeErrorPlot(ts, ener, dener, 4, "MF", "Energy", "T", "E", "MF/energy.png")
    # MakeErrorPlot(ts, magn, dmagn, 4, "MF", "Magnetization", "T", "M", "MF/magnetization.png")



    # save mean field sample
    # save("data/MeanField.jld", "ts", ts, "idx", idx, "nsamples", nsamples, "magn_sample", magn_sample, "ener_sample", ener_sample)

    

end