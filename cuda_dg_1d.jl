# Remove it after first run to avoid recompilation
#include("header.jl") 

# Use the target test header file
#= include("test/linear_scalar_advection_1d.jl") =#
include("test/compressible_euler_1d.jl")

# Kernel configurators 
#################################################################################

# CUDA kernel configurator for 1D array computing
function configurator_1d(kernel::CUDA.HostKernel, array::CuArray{Float32,1})
    config = launch_configuration(kernel.fun)

    threads = min(length(array), config.threads)
    blocks = cld(length(array), threads)

    return (threads=threads, blocks=blocks)
end

# CUDA kernel configurator for 2D array computing
function configurator_2d(kernel::CUDA.HostKernel, array::CuArray{Float32,2})
    config = launch_configuration(kernel.fun)

    threads = Tuple(fill(Int(floor((min(maximum(size(array)), config.threads))^(1 / 2))), 2))
    blocks = map(cld, size(array), threads)

    return (threads=threads, blocks=blocks)
end

# CUDA kernel configurator for 3D array computing
function configurator_3d(kernel::CUDA.HostKernel, array::CuArray{Float32,3})
    config = launch_configuration(kernel.fun)

    threads = Tuple(fill(Int(floor((min(maximum(size(array)), config.threads))^(1 / 3))), 3))
    blocks = map(cld, size(array), threads)

    return (threads=threads, blocks=blocks)
end

# Helper functions
#################################################################################

# Rewrite `get_node_vars()` as a helper function
@inline function get_nodes_vars(u, equations, indices...)

    SVector(ntuple(@inline(v -> u[v, indices...]), Val(nvariables(equations))))
end

# Rewrite `get_surface_node_vars()` as a helper function
@inline function get_surface_node_vars(u, equations, indices...)

    u_ll = SVector(ntuple(@inline(v -> u[1, v, indices...]), Val(nvariables(equations))))
    u_rr = SVector(ntuple(@inline(v -> u[2, v, indices...]), Val(nvariables(equations))))

    return u_ll, u_rr
end

# Rewrite `get_node_coords()` as a helper function
@inline function get_node_coords(x, equations, indices...)

    SVector(ntuple(@inline(idx -> x[idx, indices...]), Val(ndims(equations))))
end

# CUDA kernels 
#################################################################################

# Copy data to GPU (run as Float32)
function copy_to_gpu!(du, u)
    du = CUDA.zeros(size(du))
    u = CuArray{Float32}(u)

    return (du, u)
end

# Copy data to CPU (back to Float64)
function copy_to_cpu!(du, u)
    du = Array{Float64}(du)
    u = Array{Float64}(u)

    return (du, u)
end

# CUDA kernel for calculating fluxes along normal direction 1 
function flux_kernel!(flux_arr, u, equations::AbstractEquations{1}, flux::Function)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    k = (blockIdx().z - 1) * blockDim().z + threadIdx().z

    if (i <= size(u, 1) && j <= size(u, 2) && k <= size(u, 3))
        u_node = get_nodes_vars(u, equations, j, k)
        @inbounds flux_arr[i, j, k] = flux(u_node, 1, equations)[i]
    end

    return nothing
end

# new
function flux_kernel_new!(flux_arr, u, equations::AbstractEquations{1}, flux::Function)
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    k = (blockIdx().z - 1) * blockDim().z + threadIdx().z

    if (j <= size(u, 2) && k <= size(u, 3))
        u_node = get_nodes_vars(u, equations, j, k)
        f = flux(u_node, 1, equations)
        @inbounds begin
            for ii in axes(u, 1)
                flux_arr[ii, j, k] = f[ii]
            end
        end
    end

    return nothing
end

# CUDA kernel for calculating weak form
function weak_form_kernel!(du, derivative_dhat, flux_arr)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    k = (blockIdx().z - 1) * blockDim().z + threadIdx().z

    if (i <= size(du, 1) && j <= size(du, 2) && k <= size(du, 3))
        @inbounds begin
            for ii in axes(du, 2)
                du[i, j, k] += derivative_dhat[j, ii] * flux_arr[i, ii, k]
            end
        end
    end

    return nothing
end

# Calculate volume integral
function cuda_volume_integral!(du, u, mesh::TreeMesh{1},
    nonconservative_terms, equations,
    volume_integral::VolumeIntegralWeakForm, dg::DGSEM)

    derivative_dhat = CuArray{Float32}(dg.basis.derivative_dhat)
    flux_arr = similar(u)

    flux_kernel = @cuda launch = false flux_kernel!(flux_arr, u, equations, flux)
    flux_kernel(flux_arr, u, equations, flux; configurator_3d(flux_kernel, flux_arr)...)

    weak_form_kernel = @cuda launch = false weak_form_kernel!(du, derivative_dhat, flux_arr)
    weak_form_kernel(du, derivative_dhat, flux_arr; configurator_3d(weak_form_kernel, du)...)

    return nothing
end

# CUDA kernel for prolonging two interfaces
function prolong_interfaces_kernel!(interfaces_u, u)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    k = (blockIdx().z - 1) * blockDim().z + threadIdx().z

    if (i <= 2 && j <= size(u, 1) && k <= size(u, 3))
        @inbounds interfaces_u[i, j, k] = u[j, (2-i)*size(u, 2)+(i-1)*1, (2-i)*k+(i-1)*(k%size(u, 3)+1)]
    end

    return nothing
end

# Prolong solution to interfaces
function cuda_prolong2interfaces!(u, mesh::TreeMesh{1}, cache)

    interfaces_u = CuArray{Float32}(cache.interfaces.u)

    prolong_interfaces_kernel = @cuda launch = false prolong_interfaces_kernel!(interfaces_u, u)
    prolong_interfaces_kernel(interfaces_u, u; configurator_3d(prolong_interfaces_kernel, interfaces_u)...)

    cache.interfaces.u = interfaces_u  # Automatically copy back to CPU

    return nothing
end

# CUDA kernel for calculating surface fluxes 
function surface_flux_kernel!(surface_flux_arr, interfaces_u, equations::AbstractEquations{1}, surface_flux::FluxLaxFriedrichs) # ::Any?
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    k = (blockIdx().z - 1) * blockDim().z + threadIdx().z

    if (i == 1 && j <= size(interfaces_u, 2) && k <= size(interfaces_u, 3))
        u_ll, u_rr = get_surface_node_vars(interfaces_u, equations, k)
        @inbounds surface_flux_arr[i, j, k] = surface_flux(u_ll, u_rr, 1, equations)[j]
    end

    return nothing
end

# new
function surface_flux_kernel_new!(surface_flux_arr, interfaces_u, equations::AbstractEquations{1}, surface_flux::FluxLaxFriedrichs) # ::Any?
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    k = (blockIdx().z - 1) * blockDim().z + threadIdx().z

    if (j <= size(interfaces_u, 2) && k <= size(interfaces_u, 3))
        u_ll, u_rr = get_surface_node_vars(interfaces_u, equations, k)
        @inbounds surface_flux_arr[1, j, k] = surface_flux(u_ll, u_rr, 1, equations)[j]
    end

    return nothing
end

# CUDA kernel for setting interface fluxes
function interface_flux_kernel!(surface_flux_values, surface_flux_arr)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    k = (blockIdx().z - 1) * blockDim().z + threadIdx().z

    if (i <= size(surface_flux_values, 1) && j <= 2 && k <= size(surface_flux_values, 3))
        @inbounds surface_flux_values[i, j, k] = surface_flux_arr[1, i,
            (j-1)*k+(2-j)*((k-1)%size(surface_flux_values, 3)+iszero(k - 1)*size(surface_flux_values, 3))]
    end

    return nothing
end

# Calculate interface fluxes
function cuda_interface_flux!(mesh::TreeMesh{1}, nonconservative_terms::False,
    equations, dg::DGSEM, cache)

    surface_flux = dg.surface_integral.surface_flux
    interfaces_u = CuArray{Float32}(cache.interfaces.u)
    surface_flux_values = CuArray{Float32}(cache.elements.surface_flux_values)
    surface_flux_arr = CuArray{Float32}(undef, (1, size(interfaces_u, 2), size(interfaces_u, 3)))

    surface_flux_kernel = @cuda launch = false surface_flux_kernel!(surface_flux_arr, interfaces_u, equations, surface_flux)
    surface_flux_kernel(surface_flux_arr, interfaces_u, equations, surface_flux; configurator_3d(surface_flux_kernel, surface_flux_arr)...)

    interface_flux_kernel = @cuda launch = false interface_flux_kernel!(surface_flux_values, surface_flux_arr)
    interface_flux_kernel(surface_flux_values, surface_flux_arr; configurator_3d(interface_flux_kernel, surface_flux_values)...)

    cache.elements.surface_flux_values = surface_flux_values # Automatically copy back to CPU

    return nothing
end

# Prolong solution to boundaries
# Calculate boundary fluxes

# CUDA kernel for calculating surface integrals
function surface_integral_kernel!(du, factor_arr, surface_flux_values)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    k = (blockIdx().z - 1) * blockDim().z + threadIdx().z

    if (i <= size(du, 1) && (j == 1 || j == size(du, 2)) && k <= size(du, 3))
        @inbounds du[i, j, k] = du[i, j, k] + (-1)^isone(j) *
                                              factor_arr[isone(j)*1+(1-isone(j))*2] *
                                              surface_flux_values[i, isone(j)*1+(1-isone(j))*2, k]
    end

    return nothing
end

# Calculate surface integrals
function cuda_surface_integral!(du, mesh::TreeMesh{1}, dg::DGSEM, cache)

    factor_arr = CuArray{Float32}([dg.basis.boundary_interpolation[1, 1], dg.basis.boundary_interpolation[end, 2]])
    surface_flux_values = CuArray{Float32}(cache.elements.surface_flux_values)

    surface_integral_kernel = @cuda launch = false surface_integral_kernel!(du, factor_arr, surface_flux_values)
    surface_integral_kernel(du, factor_arr, surface_flux_values; configurator_3d(surface_integral_kernel, du)...)

    return nothing
end

# CUDA kernel for applying inverse Jacobian 
function jacobian_kernel!(du, inverse_jacobian)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    k = (blockIdx().z - 1) * blockDim().z + threadIdx().z

    if (i <= size(du, 1) && j <= size(du, 2) && k <= size(du, 3))
        @inbounds du[i, j, k] *= -inverse_jacobian[k]
    end

    return nothing
end

# Apply Jacobian from mapping to reference element
function cuda_jacobian!(du, mesh::TreeMesh{1}, cache)

    inverse_jacobian = CuArray{Float32}(cache.elements.inverse_jacobian)

    jacobian_kernel = @cuda launch = false jacobian_kernel!(du, inverse_jacobian)
    jacobian_kernel(du, inverse_jacobian; configurator_3d(jacobian_kernel, du)...)

    return nothing
end

# CUDA kernel for calculating source terms
function source_terms_kernel!(du, u, node_coordinates, t, equations::AbstractEquations{1}, source_terms::Function)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    k = (blockIdx().z - 1) * blockDim().z + threadIdx().z

    if (i <= size(du, 1) && j <= size(du, 2) && k <= size(du, 3))
        u_local = get_nodes_vars(u, equations, j, k)
        x_local = get_node_coords(node_coordinates, equations, j, k)
        @inbounds du[i, j, k] += source_terms(u_local, x_local, t, equations)[i]
    end

    return nothing
end

# Calculate source terms               
function cuda_sources!(du, u, t, source_terms::Nothing,
    equations::AbstractEquations{1}, cache)

    return nothing
end

# Calculate source terms 
function cuda_sources!(du, u, t, source_terms,
    equations::AbstractEquations{1}, cache)

    node_coordinates = CuArray{Float32}(cache.elements.node_coordinates)

    source_terms_kernel = @cuda launch = false source_terms_kernel!(du, u, node_coordinates, t, equations, source_terms)
    source_terms_kernel(du, u, node_coordinates, t, equations, source_terms; configurator_3d(source_terms_kernel, du)...)

    return nothing
end

# Inside `rhs!()` raw implementation
#################################################################################
du, u = copy_to_gpu!(du, u)

#= cuda_volume_integral!(
    du, u, mesh,
    have_nonconservative_terms(equations), equations,
    solver.volume_integral, solver)

cuda_prolong2interfaces!(u, mesh, cache)

cuda_interface_flux!(
    mesh, have_nonconservative_terms(equations),
    equations, solver, cache,)

cuda_surface_integral!(du, mesh, solver, cache)

cuda_jacobian!(du, mesh, cache)

cuda_sources!(du, u, t,
    source_terms, equations, cache)

du, u = copy_to_cpu!(du, u) =#

derivative_dhat = CuArray{Float32}(solver.basis.derivative_dhat)
flux_arr = similar(u)

@benchmark begin
    flux_kernel = @cuda launch = false flux_kernel!(flux_arr, u, equations, flux)
    flux_kernel(flux_arr, u, equations, flux; configurator_3d(flux_kernel, flux_arr)...)
end

size_arr = CuArray{Float32}(undef, (size(flux_arr, 2), size(flux_arr, 3)))
@benchmark begin
    flux_kernel_new = @cuda launch = false flux_kernel_new!(flux_arr, u, equations, flux)
    flux_kernel_new(flux_arr, u, equations, flux; configurator_2d(flux_kernel_new, size_arr)...)
end

surface_flux = solver.surface_integral.surface_flux
interfaces_u = CuArray{Float32}(cache.interfaces.u)
surface_flux_values = CuArray{Float32}(cache.elements.surface_flux_values)
surface_flux_arr = CuArray{Float32}(undef, (1, size(interfaces_u, 2), size(interfaces_u, 3)))

@benchmark begin
    surface_flux_kernel = @cuda launch = false surface_flux_kernel!(surface_flux_arr, interfaces_u, equations, surface_flux)
    surface_flux_kernel(surface_flux_arr, interfaces_u, equations, surface_flux; configurator_3d(surface_flux_kernel, surface_flux_arr)...)
end

size_arr = CuArray{Float32}(undef, (size(interfaces_u, 2), size(interfaces_u, 3)))
@benchmark begin
    surface_flux_kernel_new = @cuda launch = false surface_flux_kernel_new!(surface_flux_arr, interfaces_u, equations, surface_flux)
    surface_flux_kernel_new(surface_flux_arr, interfaces_u, equations, surface_flux; configurator_2d(surface_flux_kernel_new, size_arr)...)
end

# For tests
#################################################################################
#= reset_du!(du, solver, cache)

calc_volume_integral!(
    du, u, mesh,
    have_nonconservative_terms(equations), equations,
    solver.volume_integral, solver, cache)

prolong2interfaces!(
    cache, u, mesh, equations, solver.surface_integral, solver)

calc_interface_flux!(
    cache.elements.surface_flux_values, mesh,
    have_nonconservative_terms(equations), equations,
    solver.surface_integral, solver, cache)

calc_surface_integral!(
    du, u, mesh, equations, solver.surface_integral, solver, cache)

apply_jacobian!(
    du, mesh, equations, solver, cache)

calc_sources!(du, u, t, source_terms, equations, solver, cache) =#


#= function rhs!(du, u, t,
    mesh::TreeMesh{1}, equations,
    initial_condition, boundary_conditions, source_terms::Source,
    dg::DG, cache) where {Source}

    reset_du!(du, solver, cache)

    calc_volume_integral!(
        du, u, mesh,
        have_nonconservative_terms(equations), equations,
        solver.volume_integral, solver, cache)

    prolong2interfaces!(
        cache, u, mesh, equations, solver.surface_integral, solver)

    calc_interface_flux!(
        cache.elements.surface_flux_values, mesh,
        have_nonconservative_terms(equations), equations,
        solver.surface_integral, solver, cache)

    calc_surface_integral!(
        du, u, mesh, equations, solver.surface_integral, solver, cache)

    apply_jacobian!(
        du, mesh, equations, solver, cache)

    calc_sources!(du, u, t, source_terms, equations, dg, cache)

    return nothing
end

function semidiscretize(semi::AbstractSemidiscretization, tspan)
    u0_ode = compute_coefficients(first(tspan), semi)

    iip = true
    specialize = SciMLBase.FullSpecialize
    return ODEProblem{iip,specialize}(rhs!, u0_ode, tspan, semi)
end =#


