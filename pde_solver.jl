# dande | sesshomaru-engine
# implicit finite difference for bs-pde.
# purely vectorized, no slow loops.

using LinearAlgebra, SparseArrays

function solve_bs_implicit(S_max::Float64, K::Float64, T::Float64, r::Float64, sigma::Float64, M::Int, N::Int)
    dt, dS = T/N, S_max/M
    S = range(0, S_max, length=M+1)
    
    # tridiagonal matrix setup for implicit scheme
    j = 1:M-1
    alpha = @. 0.5 * dt * (sigma^2 * j^2 - r * j)
    beta  = @. 1.0 + dt * (sigma^2 * j^2 + r)
    gamma = @. 0.5 * dt * (sigma^2 * j^2 + r * j)
    
    A = spdiagm(-1 => -alpha[2:end], 0 => beta, 1 => -gamma[1:end-1])
    
    # terminal condition (european call)
    V = max.(S[2:end-1] .- K, 0.0)
    
    # time stepping
    for i in 1:N
        V[end] += gamma[end] * (S_max - K * exp(-r * i * dt)) # boundary
        V = A \ V
    end
    return V[Int(floor(M/2))] # return at-the-money price
end

@time price = solve_bs_implicit(200.0, 100.0, 1.0, 0.05, 0.2, 1000, 1000)
println("ATM Price: ", round(price, digits=4))
