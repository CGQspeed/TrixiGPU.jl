#= include("../cuda_dg_2d.jl") =#

equations = CompressibleEulerEquations2D(1.4)

function initial_condition_isentropic_vortex(x, t, equations::CompressibleEulerEquations2D)

    inicenter = SVector(0.0, 0.0)

    iniamplitude = 5.0

    rho = 1.0
    v1 = 1.0
    v2 = 1.0
    vel = SVector(v1, v2)
    p = 25.0
    rt = p / rho
    t_loc = 0.0
    cent = inicenter + vel * t_loc

    cent = x - cent

    cent = SVector(-cent[2], cent[1])
    r2 = cent[1]^2 + cent[2]^2
    du = iniamplitude / (2 * π) * exp(0.5 * (1 - r2))
    dtemp = -(equations.gamma - 1) / (2 * equations.gamma * rt) * du^2
    rho = rho * (1 + dtemp)^(1 / (equations.gamma - 1))
    vel = vel + du * cent
    v1, v2 = vel
    p = p * (1 + dtemp)^(equations.gamma / (equations.gamma - 1))
    prim = SVector(rho, v1, v2, p)
    return prim2cons(prim, equations)
end

initial_condition = initial_condition_isentropic_vortex
solver = DGSEM(polydeg = 3, surface_flux = flux_lax_friedrichs)

coordinates_min = (-10.0, -10.0)
coordinates_max = (10.0, 10.0)
mesh = TreeMesh(
    coordinates_min,
    coordinates_max,
    initial_refinement_level = 4,
    n_cells_max = 10_000,
)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)

tspan = (0.0, 20.0)

# Run on CPU
#################################################################################
ode_cpu = semidiscretize_cpu(semi, tspan)

sol_cpu = OrdinaryDiffEq.solve(
    ode_cpu,
    RDPK3SpFSAL49();
    abstol = 1.0e-6,
    reltol = 1.0e-6,
    ode_default_options()...,
)

# Run on GPU
#################################################################################
ode_gpu = semidiscretize_gpu(semi, tspan)

sol_gpu = OrdinaryDiffEq.solve(
    ode_gpu,
    RDPK3SpFSAL49();
    abstol = 1.0e-6,
    reltol = 1.0e-6,
    ode_default_options()...,
)

# Compare results
################################################################################
extrema(sol_cpu.u[end] - sol_gpu.u[end])
