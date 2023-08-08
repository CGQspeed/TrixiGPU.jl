# The header part for testing boundaries in 1D
equations = HyperbolicDiffusionEquations1D()

initial_condition = initial_condition_poisson_nonperiodic

boundary_conditions = boundary_condition_poisson_nonperiodic

solver = DGSEM(polydeg=4, surface_flux=flux_lax_friedrichs)

coordinates_min = 0.0
coordinates_max = 1.0
mesh = TreeMesh(coordinates_min, coordinates_max,
    initial_refinement_level=3,
    n_cells_max=30_000,
    periodicity=false)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
    boundary_conditions=boundary_conditions,
    source_terms=source_terms_poisson_nonperiodic)

@unpack mesh, equations, initial_condition, boundary_conditions, source_terms, solver, cache = semi

t = 0.0
tspan = (0.0, 5.0)

ode = semidiscretize(semi, tspan)
u_ode = copy(ode.u0)
du_ode = similar(u_ode)
u = wrap_array(u_ode, mesh, equations, solver, cache)
du = wrap_array(du_ode, mesh, equations, solver, cache)