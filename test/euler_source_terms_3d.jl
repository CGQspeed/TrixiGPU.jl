# The header part for testing source terms in 3D 
equations = CompressibleEulerEquations3D(1.4f0)

initial_condition = initial_condition_convergence_test

solver = DGSEM(polydeg=3, surface_flux=flux_lax_friedrichs,
    volume_integral=VolumeIntegralWeakForm())

coordinates_min = (0.0, 0.0, 0.0)
coordinates_max = (2.0, 2.0, 2.0)
mesh = TreeMesh(coordinates_min, coordinates_max,
    initial_refinement_level=2,
    n_cells_max=10_000)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
    source_terms=source_terms_convergence_test)

@unpack mesh, equations, initial_condition, boundary_conditions, source_terms, solver, cache = semi

t = 0.0
l = nvariables(equations) * nnodes(solver)^ndims(mesh) * nelements(solver, cache)
du_ode = ones(Float64, l)
u_ode = ones(Float64, l)
du = wrap_array(du_ode, mesh, equations, solver, cache)
u = wrap_array(u_ode, mesh, equations, solver, cache)