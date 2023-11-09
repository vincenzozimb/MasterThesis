using QuadGK
# using ForwardDiff
    

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

