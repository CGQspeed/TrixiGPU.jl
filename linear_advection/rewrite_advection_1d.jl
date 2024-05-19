### This rewriting of discretization of 1D advection equation are based on the rewrite demo
### and CUDA.jl is not yet introduced.

using Trixi, LinearAlgebra, OrdinaryDiffEq

coordinates_min = -1.0
coordinates_max = 1.0
n_elements = 16
dx = (coordinates_max - coordinates_min) / n_elements

polydeg = 3
basis = LobattoLegendreBasis(polydeg)
nodes = basis.nodes
D = basis.derivative_matrix
M = diagm(basis.weights)
B = diagm([-1; zeros(polydeg - 1); 1])

# --------------------------------------------------------------------- Rewrite Done----------
# Create column vector x_{l}
x_l = collect(0:n_elements-1) .* dx + fill(dx / 2 - 1, n_elements)

# Get matrix x from column vectors x_{l} and ξ
x =
    fill(1.0, polydeg + 1) * transpose(x_l) +
    nodes * transpose(fill(1.0, n_elements)) .* (dx / 2)
# --------------------------------------------------------------------------------------------

initial_condition_sine_wave(x) = 1.0 + 0.5 * sin(pi * x)
u0 = initial_condition_sine_wave.(x)
surface_flux = flux_lax_friedrichs

# --------------------------------------------------------------------- Rewrite Done----------
function rhs!(du, u, x, t)
    du .= zero(eltype(du))
    λ = copy(du)    # Rename flux_numerical as λ for simplicity
    equations = LinearScalarAdvectionEquation1D(1.0)

    # Caculate u_{0} and u_{N}
    u_0 =
        transpose([1.0; fill(0.0, polydeg)]) *
        u *
        Matrix{Float64}(I, n_elements, n_elements)[:, [2:n_elements; 1]]
    u_N = transpose([fill(0.0, polydeg); 1.0]) * u

    # Caculate flux vector γ 
    γ = surface_flux.(u_N, u_0, 1, equations)   # Take γ as a row vector

    # Get flux matrix λ from γ and γ'
    λ = [
        γ * Matrix{Float64}(I, n_elements, n_elements)[:, [n_elements; 1:n_elements-1]]
        zeros(polydeg - 1, n_elements)
        γ
    ]

    # Reset matrix du
    du = (-(M \ B) * λ + (M \ transpose(D)) * M * u) .* (2 / dx)

    return nothing
end
# --------------------------------------------------------------------------------------------

tspan = (0.0, 2.0)
ode = ODEProblem(rhs!, u0, tspan)
sol = solve(ode, RDPK3SpFSAL49(), abstol = 1.0e-6, reltol = 1.0e-6, save_everystep = false)
