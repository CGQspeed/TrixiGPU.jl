# The header part for testing basic kernels in 3D
advection_velocity = (0.2, -0.7, 0.5)
equations = LinearScalarAdvectionEquation3D(advection_velocity)

solver = DGSEM(polydeg=3, surface_flux=flux_lax_friedrichs)

coordinates_min = (-1.0, -1.0, -1.0)
coordinates_max = (1.0, 1.0, 1.0)

mesh = TreeMesh(coordinates_min, coordinates_max, initial_refinement_level=3, n_cells_max=30_000)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_convergence_test, solver)

@unpack mesh, equations, initial_condition, boundary_conditions, source_terms, solver, cache = semi

t = 0.0

ode = semidiscretize(semi, tspan)
u_ode = copy(ode.u0)
du_ode = similar(u_ode)
u = wrap_array(u_ode, mesh, equations, solver, cache)
du = wrap_array(du_ode, mesh, equations, solver, cache)
