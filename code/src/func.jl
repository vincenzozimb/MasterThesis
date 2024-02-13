## ---------- GENERIC FUNCTIONS AND VARIABLES ---------- ##

using QuadGK
using Plots
using LaTeXStrings


## global variables
Tc = 2.0 / (log(1.0+sqrt(2.0))) # Tc ≈ 2.2691853 / J (remember to multiply by the coupling constant)
βc = 1 / Tc
q = 4 # coordination number


## specific heat capacity with finite differences
function SpecificHeat(F, T)
    delta_T = diff(T)
    dF_dT = diff(F) ./ delta_T
    d2F_dT2 = diff(dF_dT) ./ delta_T[1:end-1]
    res = -T[1:end-2] .* d2F_dT2
    append!(res,[NaN, NaN])    
    return res
end


## moving average filter
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


## square lattice Ising exact free energy
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


## square lattice Ising exact magnetization
function ising_magnetization(β::Real, J::Real=1.0)
    # βc = log(1+sqrt(2))/2
    β > βc && return (1 - sinh(2 * β*J)^(-4))^(1 / 8)
    return 0.0
end


# add image folder to path
images_path = pwd() * "/images"
if !isdir(images_path)
    mkdir(images_path)
end


## plot function, with critical temperature
function MakePlot(x1, y1, x2, y2, Tc, lab1::String, lab2::String, title::String, xlab::String, ylab::String, saveas)
    scatter(x1, y1, ms=2, label=lab1)
    vline!([Tc], line=:red, label=L"T_c")   
    !isnan(x2[1]) ? plot!(x2, y2, label=lab2) : nothing
    # !isnan(x2[1]) ? scatter!(x2, y2, ms=2, label=lab2) : nothing
    title!(title)
    xlabel!(xlab)
    ylabel!(ylab)
    !ismissing(saveas) ? savefig(joinpath(images_path, saveas)) : display(plot!())
end


## plot with errorbar, with critical temperature
function MakeErrorPlot(x, y, dy, Tc, lab::String, title::String, xlab::String, ylab::String, saveas)
    scatter(x, y, ms=2, label=lab)
    yerror!(x, y, yerr=dy, msc=:blue)
    vline!([Tc], line=:red, label=L"T_c")
    title!(title)
    xlabel!(xlab)
    ylabel!(ylab)
    !ismissing(saveas) ? savefig(joinpath(images_path, saveas)) : display(plot!()) 
end


## Newton-Raphson method
function NewtonRaphson(f::Function, x0::Number, fprime::Function, args::Tuple=(); 
                tol::AbstractFloat=1e-8, maxiter::Integer=50, eps::AbstractFloat=1e-10)
    for _ in 1:maxiter
        yprime = fprime(x0, args...)
        if abs(yprime) < eps
            warn("First derivative is zero")
            return x0
        end
        y = f(x0, args...)
        x1 = x0 - y/yprime
        if abs(x1-x0) < tol
            return x1
        end
        x0 = x1
    end
    error("Max iteration exceeded")
end


## function for self consistent equation in mean field
function mf(L::Real, β::Real, J::Real=1.0, h::Real=0.0)
    return L - tanh(β * (h + q*J*L) )
end


## derivative of the function for self consistent equation in mean field
function mfder(L::Real, β::Real, J::Real=1.0, h::Real=0.0)
    return 1.0 - β*J*q / (cosh( β * (h + q*J*L) )^2)
end


## sample configuration from the mean field distribution
function sampleMF(Nspins::Integer, p::Real)
    n = Int(sqrt(Nspins))
    spins = rand(Bernoulli(p), Nspins)
    spins = (2 * spins) .- 1
    spins = reshape(spins, n, n) 
    return spins
end


## Ising energy functions
function IsingEnergy(spins::Array{Int, 2}; J::Float64=1.0, h::Float64=0.0)
    Lx, Ly = size(spins)
    energy = 0
    for i in 1:Lx
        for j in 1:Ly
            s = spins[i, j]
            nb_sum = spins[mod1(i+1, Lx), j] +
                        spins[mod1(i-1, Lx), j] +
                        spins[i, mod1(j+1, Ly)] +
                        spins[i, mod1(j-1, Ly)]

            energy += -J * s * nb_sum - h * s
        end
    end
    energy /= 2  # divide by 2 to avoid double counting
    return energy
end


## Mean Field magnetization
function MFmagnetization(ts, J=1.0, h=0.0)
    # (NaN at the MF critical point)
    Tcmf = 4.0
    x0 = 1.0
    L0 = []
    for t in ts
        t != Tcmf ? y = NewtonRaphson(mf, x0, mfder, (1.0/t, J, h)) : y = NaN 
        push!(L0,y)
    end
    return L0
end
