## ---------- GENERIC FUNCTIONS AND VARIABLES ---------- ##

using QuadGK
using Plots
using LaTeXStrings


## global variables
Tc = 2.0 / (log(1.0+sqrt(2.0))) # Tc ≈ 2.2691853 / J (remember to multiply by the coupling constant)
βc = 1 / Tc


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
function MakePlot(x1, y1, x2, y2, Tc, lab1::String, lab2::String, title::String, xlab::String, ylab::String, saveas::String)
    scatter(x1, y1, ms=2, label=lab1)
    vline!([Tc], line=:red, label=L"T_c")   
    !isnan(x2[1]) ? plot!(x2, y2, label=lab2) : nothing
    title!(title)
    xlabel!(xlab)
    ylabel!(ylab)
    savefig(joinpath(images_path, saveas))  
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
