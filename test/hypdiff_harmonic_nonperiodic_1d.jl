# The header part for testing boundaries in 1D
equations = HyperbolicDiffusionEquations1D(nu=1.25)

function initial_condition_harmonic_nonperiodic(x, t, equations::HyperbolicDiffusionEquations1D)
    if t == 0.0
        phi = 5.0
        q1 = 0.0
    else
        A = 3
        B = exp(1)
        phi = A + B * x[1]
        q1 = B
    end
    return SVector(phi, q1)
end

initial_condition = initial_condition_harmonic_nonperiodic

boundary_conditions = BoundaryConditionDirichlet(initial_condition)

solver = DGSEM(polydeg=3, surface_flux=flux_lax_friedrichs)

coordinates_min = -1.0
coordinates_max = 2.0
mesh = TreeMesh(coordinates_min, coordinates_max,
    initial_refinement_level=2,
    n_cells_max=30_000,
    periodicity=false)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
    boundary_conditions=boundary_conditions,
    source_terms=source_terms_harmonic)

@unpack mesh, equations, initial_condition, boundary_conditions, source_terms, solver, cache = semi

t = 0.0

ode = semidiscretize(semi, tspan)
u_ode = copy(ode.u0)
du_ode = similar(u_ode)
u = wrap_array(u_ode, mesh, equations, solver, cache)
du = wrap_array(du_ode, mesh, equations, solver, cache)