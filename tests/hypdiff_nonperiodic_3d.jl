# The header part for testing boundaries in 3D
equations = HyperbolicDiffusionEquations3D()

initial_condition = initial_condition_poisson_nonperiodic
boundary_conditions = (x_neg=boundary_condition_poisson_nonperiodic,
    x_pos=boundary_condition_poisson_nonperiodic,
    y_neg=boundary_condition_periodic,
    y_pos=boundary_condition_periodic,
    z_neg=boundary_condition_periodic,
    z_pos=boundary_condition_periodic)

solver = DGSEM(polydeg=4, surface_flux=flux_lax_friedrichs)

coordinates_min = (0.0f0, 0.0f0, 0.0f0)
coordinates_max = (1.0f0, 1.0f0, 1.0f0)
mesh = TreeMesh(coordinates_min, coordinates_max,
    initial_refinement_level=2,
    n_cells_max=30_000,
    periodicity=(false, true, true))

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
    source_terms=source_terms_poisson_nonperiodic,
    boundary_conditions=boundary_conditions)

@unpack mesh, equations, initial_condition, boundary_conditions, source_terms, solver, cache = semi

t = 0.0f0
tspan = (0.0f0, 5.0f0)

ode = semidiscretize(semi, tspan)
u_ode = copy(ode.u0)
du_ode = similar(u_ode)
u = wrap_array(u_ode, mesh, equations, solver, cache)
du = wrap_array(du_ode, mesh, equations, solver, cache)