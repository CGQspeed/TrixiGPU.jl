#= include("../cuda_dg_2d.jl") =#

advection_velocity = (0.2, -0.7)
equations = LinearScalarAdvectionEquation2D(advection_velocity)

initial_condition = initial_condition_convergence_test
solver = DGSEM(polydeg=3, surface_flux=flux_lax_friedrichs)

coordinates_min = (-1.0, -1.0)
coordinates_max = (1.0, 1.0)
refinement_patches = (
    (type="box", coordinates_min=(0.0, -1.0), coordinates_max=(1.0, 1.0)),
)
mesh = TreeMesh(coordinates_min, coordinates_max,
    initial_refinement_level=2,
    refinement_patches=refinement_patches,
    n_cells_max=10_000,)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)

tspan = (0.0, 1.0)

# Run on CPU
#################################################################################
ode_cpu = semidiscretize_cpu(semi, tspan)

sol_cpu = OrdinaryDiffEq.solve(ode_cpu, RDPK3SpFSAL49();
    abstol=1.0e-6, reltol=1.0e-6, ode_default_options()...)

#= u0_ode_cpu = copy(ode_cpu.u0)
du_ode_cpu = similar(u0_ode_cpu)
Trixi.rhs!(du_ode_cpu, u0_ode_cpu, semi, 0.0) =#

# Run on GPU
#################################################################################
ode_gpu = semidiscretize_gpu(semi, tspan)

sol_gpu = OrdinaryDiffEq.solve(ode_gpu, RDPK3SpFSAL49();
    abstol=1.0e-6, reltol=1.0e-6, ode_default_options()...)

#= u0_ode_gpu = copy(ode_gpu.u0)
du_ode_gpu = similar(u0_ode_gpu)
Trixi.rhs!(du_ode_gpu, u0_ode_gpu, semi, 0.0) =#

# Compare results
################################################################################
extrema(sol_cpu.u[end] - sol_gpu.u[end])