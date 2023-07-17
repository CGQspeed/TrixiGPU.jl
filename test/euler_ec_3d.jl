# The header part for testing flux differencing in 3D
equations = CompressibleEulerEquations3D(1.4)

initial_condition = initial_condition_weak_blast_wave

volume_flux = flux_ranocha
solver = DGSEM(polydeg=3, surface_flux=flux_ranocha,
    volume_integral=VolumeIntegralFluxDifferencing(volume_flux))

coordinates_min = (-2.0, -2.0, -2.0)
coordinates_max = (2.0, 2.0, 2.0)
mesh = TreeMesh(coordinates_min, coordinates_max,
    initial_refinement_level=3,
    n_cells_max=100_000)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)

@unpack mesh, equations, initial_condition, boundary_conditions, source_terms, solver, cache = semi

t = 0.0
l = nvariables(equations) * nnodes(solver)^ndims(mesh) * nelements(solver, cache)
du_ode = ones(Float64, l)
u_ode = ones(Float64, l)
du = wrap_array(du_ode, mesh, equations, solver, cache)
u = wrap_array(u_ode, mesh, equations, solver, cache)