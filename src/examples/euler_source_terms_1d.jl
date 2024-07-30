#= include("../cuda_dg_1d.jl") =#

# Run on CPU
#################################################################################
equations = CompressibleEulerEquations1D(1.4)

initial_condition = initial_condition_convergence_test

solver = DGSEM(polydeg = 4, surface_flux = flux_lax_friedrichs)

coordinates_min = 0.0
coordinates_max = 2.0

mesh = TreeMesh(coordinates_min,
                coordinates_max,
                initial_refinement_level = 4,
                n_cells_max = 10_000)

semi = SemidiscretizationHyperbolic(mesh,
                                    equations,
                                    initial_condition,
                                    solver,
                                    source_terms = source_terms_convergence_test)

tspan = (0.0, 2.0)

ode_cpu = semidiscretize_cpu(semi, tspan)

sol_cpu = OrdinaryDiffEq.solve(ode_cpu,
                               BS3(); #= SSPRK43() =#
                               abstol = 1.0e-6,
                               reltol = 1.0e-6,
                               ode_default_options()...,)

# Run on GPU
#################################################################################
equations = CompressibleEulerEquations1D(1.4f0)

initial_condition = initial_condition_convergence_test

solver = DGSEM(polydeg = 4, surface_flux = flux_lax_friedrichs)

coordinates_min = 0.0f0
coordinates_max = 2.0f0

mesh = TreeMesh(coordinates_min,
                coordinates_max,
                initial_refinement_level = 4,
                n_cells_max = 10_000)

semi = SemidiscretizationHyperbolic(mesh,
                                    equations,
                                    initial_condition,
                                    solver,
                                    source_terms = source_terms_convergence_test)

tspan = (0.0f0, 2.0f0)

ode_gpu = semidiscretize_gpu(semi, tspan)

sol_gpu = OrdinaryDiffEq.solve(ode_gpu,
                               BS3(); #= SSPRK43() =#
                               abstol = 1.0e-6,
                               reltol = 1.0e-6,
                               ode_default_options()...,)

# Compare results
################################################################################
extrema(sol_cpu.u[end] - sol_gpu.u[end])
