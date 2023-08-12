# The header part for testing true nonconservative terms in 2D
equations = ShallowWaterEquations1D(gravity_constant=9.81, H0=3.25)

function initial_condition_discontinuous_well_balancedness(x, t, equations::ShallowWaterEquations1D)

    H = equations.H0
    v = 0.0
    b = 0.0

    if x[1] >= 0.5 && x[1] <= 0.75
        b = 2.0 + 0.5 * sin(2.0 * pi * x[1])
    end

    return prim2cons(SVector(H, v, b), equations)
end

initial_condition = initial_condition_discontinuous_well_balancedness

volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
surface_flux = (flux_fjordholm_etal, flux_nonconservative_fjordholm_etal)
solver = DGSEM(polydeg=4, surface_flux=surface_flux,
    volume_integral=VolumeIntegralFluxDifferencing(volume_flux))

coordinates_min = -1.0
coordinates_max = 1.0
mesh = TreeMesh(coordinates_min, coordinates_max,
    initial_refinement_level=3,
    n_cells_max=10_000)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)

@unpack mesh, equations, initial_condition, boundary_conditions, source_terms, solver, cache = semi

t = 0.0
tspan = (0.0, 100.0)

ode = semidiscretize(semi, tspan)
u_ode = copy(ode.u0)
du_ode = similar(u_ode)
u = wrap_array(u_ode, mesh, equations, solver, cache)
du = wrap_array(du_ode, mesh, equations, solver, cache)