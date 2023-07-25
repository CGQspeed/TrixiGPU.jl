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
ode = semidiscretize_gpu(semi, tspan)

sol = OrdinaryDiffEq.solve(ode, SSPRK43();
    abstol=1.0e-6, reltol=1.0e-6, ode_default_options()...)