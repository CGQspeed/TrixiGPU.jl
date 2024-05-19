### Introduce CUDA.jl into the implementation of discretization of 1D advection equation. 
### Note:
### 1. This raw implementation is not yet optimized;
### 2. All the numbers are converted to the type of Float32;
### 3. One custom kernel is applied and the rest operations are based on built-in CUDA functions.

using Trixi, LinearAlgebra, OrdinaryDiffEq, CUDA, Plots

# Get CPU scalar data
coordinates_min = -1.0
coordinates_max = 1.0
n_elements = 16
polydeg = 3
dx = (coordinates_max - coordinates_min) / n_elements

# Create basis with Lagrange polynomials
basis = LobattoLegendreBasis(polydeg)

# --------------------------------------------------------------------- GPU Acceleration----------    
# Get CPU array data and copy to GPU
nodes = CuArray{Float32}(basis.nodes)
D = CuArray(convert(Array{Float32}, basis.derivative_matrix))
M = CuArray(convert(Array{Float32}, diagm(basis.weights)))
B = CuArray(convert(Array{Float32}, diagm([-1; zeros(polydeg - 1); 1])))

# Create column vector x_{l} on GPU                            
x_l =
    CuArray{Float32}(collect(0:n_elements-1)) .* Float32(dx) +
    CuArray{Float32}(fill(dx / 2 - 1, n_elements))

# Get matrix x from column vectors x_{l} and ξ on GPU           
x =
    CuArray{Float32}(fill(1.0, polydeg + 1)) * transpose(x_l) +
    nodes * transpose(CuArray{Float32}(fill(1, n_elements))) .* (Float32(dx) / 2)

# Get u0 by using map function on GPU
u0 = map(x -> Float32(0.5) * sin(pi * x) + 1, x)

# Write a CUDA kernel for numerical flux
function apply_surface_flux!(γ, u1, u2)
    # Get index and stride for each thread
    index = threadIdx().x
    stride = blockDim().x

    # Choose to use Lax-Friedrichs flux and set equations
    surface_flux = flux_lax_friedrichs
    equations = LinearScalarAdvectionEquation1D(1.0)

    # Assign each thread to calculate corresponding flux
    for i ∈ index:stride:length(u1)
        @inbounds γ[i] = surface_flux(u1[i], u2[i], 1, equations)
    end

    return nothing
end

# Apply the above CUDA kernel within rhs!() function            
function rhs!(du, u, x, t)
    # Initialize u and du and copy them to GPU
    u = CuArray(convert(Array{Float32}, u))
    du = CuArray{Float32}(undef, (polydeg + 1, n_elements))

    # Caculate u_{0} and u_{N} on GPU  
    u_0 =
        transpose([1.0f0; CUDA.fill(0.0f0, polydeg)]) *
        u *
        CuArray(Matrix{Float32}(I, n_elements, n_elements)[:, [2:n_elements; 1]])
    u_N = transpose([CUDA.fill(0.0f0, polydeg); 1.0f0]) * u

    # Initialize γ and call the custom CUDA kernel
    γ = similar(u_0)
    @cuda threads = 8 apply_surface_flux!(γ, u_0, u_N)

    # Get flux matrix λ from γ and γ' on GPU  
    λ = [
        γ * CuArray(
            Matrix{Float32}(I, n_elements, n_elements)[:, [n_elements; 1:n_elements-1]],
        )
        CUDA.zeros(polydeg - 1, n_elements)
        γ
    ]

    # Reset matrix du on GPU 
    du = (-(M \ B) * λ + (M \ transpose(D)) * M * u) .* (2 / dx)

    # Copy u and du back to CPU
    u = Array(u)
    du = Array(du)

    return nothing
end

# Copy GPU data back to CPU
u0 = Array(u0)
# ------------------------------------------------------------------------------------------------

# Create and solve ODE problem
tspan = (0.0, 2.0)
ode = ODEProblem(rhs!, u0, tspan)
sol = solve(ode, RDPK3SpFSAL49(), abstol = 1.0e-6, reltol = 1.0e-6, save_everystep = false)

#= # Plot the solution
plot(vec(x), vec(sol.u[end]), label="solution at t=$(tspan[2])", legend=:topleft, lw=3) =#
