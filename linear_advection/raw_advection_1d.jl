### Implementation of discretization of 1D advection equation and then solve the 
### transformed ODE problem.

using Trixi, LinearAlgebra, OrdinaryDiffEq, Plots

# Create 1D mesh by setting minimum coordinate and maximum coordinate(physical domain),
# give the number of elements Q_{l} to split for physical domain and calculate 
# the length of each element Q_{l}
coordinates_min = -1.0
coordinates_max = 1.0
n_elements = 16
dx = (coordinates_max - coordinates_min) / n_elements

# Define the polynomial degree to approximate the solution in each Q_{l},
# the Lagrange basis funtions are chosen as the polynomial basis
polydeg = 3

# Create Gauss-Lobatto nodes with degree equals to polydeg, these nodes are 
# ξ_{0}, ..., ξ_{N} in interval [-1, 1](reference domain)
# N the is the degree of Lagrange basis function(u_{0}, ..., u_{N} and l_{0}, ..., l_{N})
# and is also the degree of Lobatto Quadrature(ξ_{0}, ..., ξ_{N})
basis = LobattoLegendreBasis(polydeg)
nodes = basis.nodes

# Create derivative matrix D, mass matrix M, and boundary matrix B
# Matrix D, M, and B are used for calculating du(t) = (du_{0}(t), ..., du_{N}(t))^{T}
# for each Q_{l} and here du is a simple representation of ∂u(x, t)/∂t 
D = basis.derivative_matrix
M = diagm(basis.weights)
B = diagm([-1; zeros(polydeg - 1); 1])

# Create a N+1 by n_elements matrix x, each element x_{ij} is the transformed 
# coordinate of ξ_{i} (i = 0, 1, ..., N) in Q_{j} (j = 1, 2, ..., n_elements)
# The transformed coordinates are in the physical domain
# --------------------------------------------------------------------- Need Rewrite--------
x = Matrix{Float64}(undef, length(nodes), n_elements)
for element = 1:n_elements
    x_l = -1 + (element - 1) * dx + dx / 2
    for i in eachindex(nodes)       # Change from length() to eachindex()
        ξ = nodes[i]
        x[i, element] = x_l + dx / 2 * ξ
    end
end
# ------------------------------------------------------------------------------------------

# Set initial condition u(x, 0) = u0(x) and take each value from matrix x into u0(x) 
# thus get the initial values for solving ODE later 
initial_condition_sine_wave(x) = 1.0 + 0.5 * sin(pi * x)
u0 = initial_condition_sine_wave.(x)

# Choose to use local Lax-Friedrichs flux to handle the different values at the same
# point at the interfaces
# The function u^{Q_{l}}(ξ, t) takes ξ_{N} = 1 has a different value compared to the 
# function u^{Q_{l+1}}(ξ, t) takes ξ_{0} = -1 at the same value t
surface_flux = flux_lax_friedrichs

# Implement a function to calculate du(t) = (du_{0}(t), ..., du_{N}(t))^{T} based on 
# the given vectors, u(t) = (u_{0}(t), ..., u_{N}(t))^{T} and u^{*}(t), for each Q_{l}
# --------------------------------------------------------------------- Need Rewrite--------
function rhs!(du, u, x, t)
    # Reset du matrix to zero
    # du is set to zero in solve() before being taken into the algo if f(rhs!) is inplace, 
    # but comment this line will cause "retcode: Unstable"
    du .= zero(eltype(du))

    # Set flux matrix to zero
    flux_numerical = copy(du)

    # Calculate interface and boundary fluxes, given function f^{*}(u) = u^{*}(u_{L}, u_{R})
    equations = LinearScalarAdvectionEquation1D(1.0)

    # In this for loop, values in flux_numerical are computed repeatedly, not both left 
    # interface and right interface are needed for each element
    for element = 2:n_elements-1

        # Get left interface of Q_{l} where l = element 
        flux_numerical[1, element] =
            surface_flux(u[end, element-1], u[1, element], 1, equations)
        flux_numerical[end, element-1] = flux_numerical[1, element]

        # Get right interface of Q_{l} where l = element
        flux_numerical[end, element] =
            surface_flux(u[end, element], u[1, element+1], 1, equations)
        flux_numerical[1, element+1] = flux_numerical[end, element]
    end

    # Calculate boundary flux, the last value of the last Q_{l} is used as u_{L} in the front 
    # boundary interface and similarly the first value of the first Q_{l} is used as u_{R}
    # in the end boundary interface
    flux_numerical[1, 1] = surface_flux(u[end, end], u[1, 1], 1, equations)
    flux_numerical[end, end] = flux_numerical[1, 1]

    # Calculate surface integrals, each element corresponds to each Q_{l}
    for element = 1:n_elements
        du[:, element] -= (M \ B) * flux_numerical[:, element]
    end

    # Calculate volume integrals, each element corresponds to each Q_{l}
    for element = 1:n_elements
        flux = u[:, element]
        du[:, element] += (M \ transpose(D)) * M * flux
    end

    # Apply Jacobian from mapping to reference element
    for element = 1:n_elements
        du[:, element] *= 2 / dx
    end

    return nothing
end
# ------------------------------------------------------------------------------------------

# Create ODE problem
tspan = (0.0, 2.0)
# No need to take in x parameter, since x does not act as the parameter outside u(x, t) 
# and the rhs!() function does not use x
ode = ODEProblem(rhs!, u0, tspan)

# Solve ODE problem
sol = solve(ode, RDPK3SpFSAL49(), abstol = 1.0e-6, reltol = 1.0e-6, save_everystep = false)

# @ODEProblem() SciMLBase.jl/src/problems/ode_problems.jl
# @solve() ODE.jl/src/common.jl
# @flux_lax_friedrichs() Trixi.jl/src/equations/numerical_fluxes.jl
