#= include("../cuda_dg_2d.jl") =#

equations = CompressibleEulerEquations2D(1.4)

initial_condition = initial_condition_convergence_test
solver = DGSEM(polydeg=3, surface_flux=flux_lax_friedrichs)

coordinates_min = (0.0, 0.0)
coordinates_max = (2.0, 2.0)
mesh = TreeMesh(coordinates_min, coordinates_max,
    initial_refinement_level=4,
    n_cells_max=10_000)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
    source_terms=source_terms_convergence_test)

tspan = (0.0, 2.0)

# Run on CPU
#################################################################################
ode_cpu = semidiscretize_cpu(semi, tspan)

sol_cpu = OrdinaryDiffEq.solve(ode_cpu, SSPRK43();
    abstol=1.0e-6, reltol=1.0e-6, ode_default_options()...)

#= u0_ode_cpu = copy(ode_cpu.u0)
du_ode_cpu = similar(u0_ode_cpu)
Trixi.rhs!(du_ode_cpu, u0_ode_cpu, semi, 0.0) =#

# Run on GPU
#################################################################################
ode_gpu = semidiscretize_gpu(semi, tspan)

sol_gpu = OrdinaryDiffEq.solve(ode_gpu, SSPRK43();
    abstol=1.0e-6, reltol=1.0e-6, ode_default_options()...)

#= u0_ode_gpu = copy(ode_gpu.u0)
du_ode_gpu = similar(u0_ode_gpu)
Trixi.rhs!(du_ode_gpu, u0_ode_gpu, semi, 0.0) =#

# Compare results
################################################################################
extrema(sol_cpu.u[end] - sol_gpu.u[end])