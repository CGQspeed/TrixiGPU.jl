# The header part for testing basic kernels in 3D
advection_velocity = (0.2f0, -0.7f0, 0.5f0)
equations = LinearScalarAdvectionEquation3D(advection_velocity)

solver = DGSEM(polydeg=3, surface_flux=flux_lax_friedrichs)

coordinates_min = (-1.0, -1.0, -1.0)
coordinates_max = (1.0, 1.0, 1.0)

mesh = TreeMesh(coordinates_min, coordinates_max, initial_refinement_level=3, n_cells_max=30_000)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_convergence_test, solver)

@unpack mesh, equations, initial_condition, boundary_conditions, source_terms, solver, cache = semi

t = 0.0
l = nvariables(equations) * nnodes(solver)^ndims(mesh) * nelements(solver, cache)
du_ode = rand(Float64, l)
u_ode = rand(Float64, l)
du = wrap_array(du_ode, mesh, equations, solver, cache)
u = wrap_array(u_ode, mesh, equations, solver, cache)
