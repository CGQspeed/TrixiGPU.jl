# Remove it after first run to avoid recompilation
#= include("header.jl") =#

# Use the target test header file
#= include("tests/advection_basic_1d.jl") =#
#= include("tests/euler_ec_1d.jl") =#
#= include("tests/euler_source_terms_1d.jl") =#
#= include("tests/hypdiff_nonperiodic_1d.jl") =#
#= include("tests/shallowwater_well_balanced_1d.jl") =#

# Kernel configurators 
#################################################################################

# CUDA kernel configurator for 1D array computing
function configurator_1d(kernel::CUDA.HostKernel, array::CuArray{<:Any,1})
    config = launch_configuration(kernel.fun)

    threads = min(length(array), config.threads)
    blocks = cld(length(array), threads)

    return (threads = threads, blocks = blocks)
end

# CUDA kernel configurator for 2D array computing
function configurator_2d(kernel::CUDA.HostKernel, array::CuArray{<:Any,2})
    config = launch_configuration(kernel.fun)

    threads =
        Tuple(fill(Int(floor((min(maximum(size(array)), config.threads))^(1 / 2))), 2))
    blocks = map(cld, size(array), threads)

    return (threads = threads, blocks = blocks)
end

# CUDA kernel configurator for 3D array computing
function configurator_3d(kernel::CUDA.HostKernel, array::CuArray{<:Any,3})
    config = launch_configuration(kernel.fun)

    threads =
        Tuple(fill(Int(floor((min(maximum(size(array)), config.threads))^(1 / 3))), 3))
    blocks = map(cld, size(array), threads)

    return (threads = threads, blocks = blocks)
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

# Helper function for checking `cache.mortars`
@inline function check_cache_mortars(cache)

    if iszero(length(cache.mortars.orientations))
        return True()
    else
        return False()
    end
end

# Helper function for stable calls to `boundary_conditions`
@generated function boundary_stable_helper(
    boundary_conditions,
    u_inner,
    orientation,
    direction,
    x,
    t,
    surface_flux,
    equations,
)

    n = length(boundary_conditions.parameters)
    quote
        @nif $n d -> d == direction d -> return boundary_conditions[d](
            u_inner,
            orientation,
            direction,
            x,
            t,
            surface_flux,
            equations,
        )
    end
end

# CUDA kernels 
#################################################################################

# Copy data to GPU (run as Float64)
function copy_to_gpu!(du, u)

    du = CUDA.zeros(Float64, size(du))
    u = CuArray{Float64}(u)

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

    j = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    k = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if (j <= size(u, 2) && k <= size(u, 3))
        u_node = get_nodes_vars(u, equations, j, k)

        flux_node = flux(u_node, 1, equations)

        @inbounds begin
            for ii in axes(u, 1)
                flux_arr[ii, j, k] = flux_node[ii]
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

# CUDA kernel for calculating volume fluxes in direction x
function volume_flux_kernel!(
    volume_flux_arr,
    u,
    equations::AbstractEquations{1},
    volume_flux::Function,
)

    j = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    k = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if (j <= size(u, 2)^2 && k <= size(u, 3))
        j1 = div(j - 1, size(u, 2)) + 1
        j2 = rem(j - 1, size(u, 2)) + 1

        u_node = get_nodes_vars(u, equations, j1, k)
        u_node1 = get_nodes_vars(u, equations, j2, k)

        volume_flux_node = volume_flux(u_node, u_node1, 1, equations)

        @inbounds begin
            for ii in axes(u, 1)
                volume_flux_arr[ii, j1, j2, k] = volume_flux_node[ii]
            end
        end
    end

    return nothing
end

# CUDA kernel for calculating symmetric and nonsymmetric fluxes in direction x
function symmetric_noncons_flux_kernel!(
    symmetric_flux_arr,
    noncons_flux_arr,
    u,
    derivative_split,
    equations::AbstractEquations{1},
    symmetric_flux::Function,
    nonconservative_flux::Function,
)

    j = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    k = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if (j <= size(u, 2)^2 && k <= size(u, 3))
        j1 = div(j - 1, size(u, 2)) + 1
        j2 = rem(j - 1, size(u, 2)) + 1

        u_node = get_nodes_vars(u, equations, j1, k)
        u_node1 = get_nodes_vars(u, equations, j2, k)

        symmetric_flux_node = symmetric_flux(u_node, u_node1, 1, equations)
        noncons_flux_node = nonconservative_flux(u_node, u_node1, 1, equations)

        @inbounds begin
            for ii in axes(u, 1)
                symmetric_flux_arr[ii, j1, j2, k] = symmetric_flux_node[ii]
                noncons_flux_arr[ii, j1, j2, k] =
                    noncons_flux_node[ii] * derivative_split[j1, j2]
            end
        end
    end

    return nothing
end

# CUDA kernel for calculating volume integrals
function volume_integral_kernel!(du, derivative_split, volume_flux_arr)

    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    k = (blockIdx().z - 1) * blockDim().z + threadIdx().z

    if (i <= size(du, 1) && j <= size(du, 2) && k <= size(du, 3))
        @inbounds begin
            for ii in axes(du, 2)
                du[i, j, k] += derivative_split[j, ii] * volume_flux_arr[i, j, ii, k]
            end
        end
    end

    return nothing
end

# CUDA kernel for calculating symmetric and nonsymmetric volume integrals
function symmetric_noncons_integral_kernel!(
    du,
    derivative_split,
    symmetric_flux_arr,
    noncons_flux_arr,
)

    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    k = (blockIdx().z - 1) * blockDim().z + threadIdx().z

    if (i <= size(du, 1) && j <= size(du, 2) && k <= size(du, 3))
        @inbounds begin
            integral_contribution = 0.0

            for ii in axes(du, 2)
                du[i, j, k] += derivative_split[j, ii] * symmetric_flux_arr[i, j, ii, k]
                integral_contribution += noncons_flux_arr[i, j, ii, k]
            end

            du[i, j, k] += 0.5 * integral_contribution
        end
    end

    return nothing
end

# Launch CUDA kernels to calculate volume integrals
function cuda_volume_integral!(
    du,
    u,
    mesh::TreeMesh{1},
    nonconservative_terms,
    equations,
    volume_integral::VolumeIntegralWeakForm,
    dg::DGSEM,
)

    derivative_dhat = CuArray{Float64}(dg.basis.derivative_dhat)
    flux_arr = similar(u)

    size_arr = CuArray{Float64}(undef, size(u, 2), size(u, 3))

    flux_kernel = @cuda launch = false flux_kernel!(flux_arr, u, equations, flux)
    flux_kernel(flux_arr, u, equations, flux; configurator_2d(flux_kernel, size_arr)...)

    weak_form_kernel = @cuda launch = false weak_form_kernel!(du, derivative_dhat, flux_arr)
    weak_form_kernel(
        du,
        derivative_dhat,
        flux_arr;
        configurator_3d(weak_form_kernel, du)...,
    )

    return nothing
end

# Launch CUDA kernels to calculate volume integrals
function cuda_volume_integral!(
    du,
    u,
    mesh::TreeMesh{1},
    nonconservative_terms::False,
    equations,
    volume_integral::VolumeIntegralFluxDifferencing,
    dg::DGSEM,
)

    volume_flux = volume_integral.volume_flux

    derivative_split = CuArray{Float64}(dg.basis.derivative_split)
    volume_flux_arr =
        CuArray{Float64}(undef, size(u, 1), size(u, 2), size(u, 2), size(u, 3))

    size_arr = CuArray{Float64}(undef, size(u, 2)^2, size(u, 3))

    volume_flux_kernel =
        @cuda launch = false volume_flux_kernel!(volume_flux_arr, u, equations, volume_flux)
    volume_flux_kernel(
        volume_flux_arr,
        u,
        equations,
        volume_flux;
        configurator_2d(volume_flux_kernel, size_arr)...,
    )

    volume_integral_kernel =
        @cuda launch = false volume_integral_kernel!(du, derivative_split, volume_flux_arr)
    volume_integral_kernel(
        du,
        derivative_split,
        volume_flux_arr;
        configurator_3d(volume_integral_kernel, du)...,
    )

    return nothing
end

# Launch CUDA kernels to calculate volume integrals
function cuda_volume_integral!(
    du,
    u,
    mesh::TreeMesh{1},
    nonconservative_terms::True,
    equations,
    volume_integral::VolumeIntegralFluxDifferencing,
    dg::DGSEM,
)

    symmetric_flux, nonconservative_flux = dg.volume_integral.volume_flux

    derivative_split = CuArray{Float64}(dg.basis.derivative_split)
    symmetric_flux_arr =
        CuArray{Float64}(undef, size(u, 1), size(u, 2), size(u, 2), size(u, 3))
    noncons_flux_arr =
        CuArray{Float64}(undef, size(u, 1), size(u, 2), size(u, 2), size(u, 3))

    size_arr = CuArray{Float64}(undef, size(u, 2)^2, size(u, 3))

    symmetric_noncons_flux_kernel = @cuda launch = false symmetric_noncons_flux_kernel!(
        symmetric_flux_arr,
        noncons_flux_arr,
        u,
        derivative_split,
        equations,
        symmetric_flux,
        nonconservative_flux,
    )
    symmetric_noncons_flux_kernel(
        symmetric_flux_arr,
        noncons_flux_arr,
        u,
        derivative_split,
        equations,
        symmetric_flux,
        nonconservative_flux;
        configurator_2d(symmetric_noncons_flux_kernel, size_arr)...,
    )

    symmetric_noncons_integral_kernel =
        @cuda launch = false symmetric_noncons_integral_kernel!(
            du,
            derivative_split,
            symmetric_flux_arr,
            noncons_flux_arr,
        )
    symmetric_noncons_integral_kernel(
        du,
        derivative_split,
        symmetric_flux_arr,
        noncons_flux_arr;
        configurator_3d(symmetric_noncons_integral_kernel, du)...,
    )

    return nothing
end

# CUDA kernel for prolonging two interfaces in direction x
function prolong_interfaces_kernel!(interfaces_u, u, neighbor_ids)

    j = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    k = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if (j <= size(interfaces_u, 2) && k <= size(interfaces_u, 3))
        left_element = neighbor_ids[1, k]
        right_element = neighbor_ids[2, k]

        @inbounds begin
            interfaces_u[1, j, k] = u[j, size(u, 2), left_element]
            interfaces_u[2, j, k] = u[j, 1, right_element]
        end
    end

    return nothing
end

# Launch CUDA kernel to prolong solution to interfaces
function cuda_prolong2interfaces!(u, mesh::TreeMesh{1}, cache)

    neighbor_ids = CuArray{Int}(cache.interfaces.neighbor_ids)
    interfaces_u = CuArray{Float64}(cache.interfaces.u)

    size_arr = CuArray{Float64}(undef, size(interfaces_u, 2), size(interfaces_u, 3))

    prolong_interfaces_kernel =
        @cuda launch = false prolong_interfaces_kernel!(interfaces_u, u, neighbor_ids)
    prolong_interfaces_kernel(
        interfaces_u,
        u,
        neighbor_ids;
        configurator_2d(prolong_interfaces_kernel, size_arr)...,
    )

    cache.interfaces.u = interfaces_u  # Automatically copy back to CPU

    return nothing
end

# CUDA kernel for calculating surface fluxes 
function surface_flux_kernel!(
    surface_flux_arr,
    interfaces_u,
    equations::AbstractEquations{1},
    surface_flux::Any,
)

    k = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    if (k <= size(surface_flux_arr, 3))
        u_ll, u_rr = get_surface_node_vars(interfaces_u, equations, k)

        surface_flux_node = surface_flux(u_ll, u_rr, 1, equations)

        @inbounds begin
            for jj in axes(surface_flux_arr, 2)
                surface_flux_arr[1, jj, k] = surface_flux_node[jj]
            end
        end
    end

    return nothing
end

# CUDA kernel for setting interface fluxes on orientation 1 
function interface_flux_kernel!(surface_flux_values, surface_flux_arr, neighbor_ids)

    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    k = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if (i <= size(surface_flux_values, 1) && k <= size(surface_flux_arr, 3))
        left_id = neighbor_ids[1, k]
        right_id = neighbor_ids[2, k]

        @inbounds begin
            surface_flux_values[i, 2, left_id] = surface_flux_arr[1, i, k]
            surface_flux_values[i, 1, right_id] = surface_flux_arr[1, i, k]
        end
    end

    return nothing
end

# Launch CUDA kernels to calculate interface fluxes
function cuda_interface_flux!(
    mesh::TreeMesh{1},
    nonconservative_terms::False,
    equations,
    dg::DGSEM,
    cache,
)

    surface_flux = dg.surface_integral.surface_flux

    neighbor_ids = CuArray{Int}(cache.interfaces.neighbor_ids)
    interfaces_u = CuArray{Float64}(cache.interfaces.u)
    surface_flux_arr = CuArray{Float64}(undef, 1, size(interfaces_u)[2:end]...)
    surface_flux_values = CuArray{Float64}(cache.elements.surface_flux_values)

    size_arr = CuArray{Float64}(undef, size(interfaces_u, 3))

    surface_flux_kernel = @cuda launch = false surface_flux_kernel!(
        surface_flux_arr,
        interfaces_u,
        equations,
        surface_flux,
    )
    surface_flux_kernel(
        surface_flux_arr,
        interfaces_u,
        equations,
        surface_flux;
        configurator_1d(surface_flux_kernel, size_arr)...,
    )

    size_arr = CuArray{Float64}(undef, size(surface_flux_values, 1), size(interfaces_u, 3))

    interface_flux_kernel = @cuda launch = false interface_flux_kernel!(
        surface_flux_values,
        surface_flux_arr,
        neighbor_ids,
    )
    interface_flux_kernel(
        surface_flux_values,
        surface_flux_arr,
        neighbor_ids;
        configurator_2d(interface_flux_kernel, size_arr)...,
    )

    cache.elements.surface_flux_values = surface_flux_values # Automatically copy back to CPU

    return nothing
end

# CUDA kernel for prolonging two boundaries in direction x
function prolong_boundaries_kernel!(boundaries_u, u, neighbor_ids, neighbor_sides)

    j = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    k = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if (j <= size(boundaries_u, 2) && k <= size(boundaries_u, 3))
        element = neighbor_ids[k]
        side = neighbor_sides[k]

        @inbounds begin
            boundaries_u[1, j, k] = u[j, size(u, 2), element] * isequal(side, 1) # set to 0 instead of NaN
            boundaries_u[2, j, k] = u[j, 1, element] * (1 - isequal(side, 1)) # set to 0 instead of NaN
        end
    end

    return nothing
end

# Assert 
function cuda_prolong2boundaries!(
    u,
    mesh::TreeMesh{1},
    boundary_condition::BoundaryConditionPeriodic,
    cache,
)

    @assert iszero(length(cache.boundaries.orientations))
end

# Launch CUDA kernel to prolong solution to boundaries
function cuda_prolong2boundaries!(
    u,
    mesh::TreeMesh{1},
    boundary_conditions::NamedTuple,
    cache,
)

    neighbor_ids = CuArray{Int}(cache.boundaries.neighbor_ids)
    neighbor_sides = CuArray{Int}(cache.boundaries.neighbor_sides)
    boundaries_u = CuArray{Float64}(cache.boundaries.u)

    size_arr = CuArray{Float64}(undef, size(boundaries_u, 2), size(boundaries_u, 3))

    prolong_boundaries_kernel = @cuda launch = false prolong_boundaries_kernel!(
        boundaries_u,
        u,
        neighbor_ids,
        neighbor_sides,
    )
    prolong_boundaries_kernel(
        boundaries_u,
        u,
        neighbor_ids,
        neighbor_sides;
        configurator_2d(prolong_boundaries_kernel, size_arr)...,
    )

    cache.boundaries.u = boundaries_u  # Automatically copy back to CPU

    return nothing
end

# CUDA kernel for getting last and first indices
function last_first_indices_kernel!(lasts, firsts, n_boundaries_per_direction)

    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    if (i <= length(n_boundaries_per_direction))
        @inbounds begin
            for ii = 1:i
                lasts[i] += n_boundaries_per_direction[ii]
            end
            firsts[i] = lasts[i] - n_boundaries_per_direction[i] + 1
        end
    end

    return nothing
end

# CUDA kernel for calculating boundary fluxes on direction 1, 2
function boundary_flux_kernel!(
    surface_flux_values,
    boundaries_u,
    node_coordinates,
    t,
    boundary_arr,
    indices_arr,
    neighbor_ids,
    neighbor_sides,
    orientations,
    boundary_conditions::NamedTuple,
    equations::AbstractEquations{1},
    surface_flux::Any,
)

    k = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    if (k <= length(boundary_arr))
        boundary = boundary_arr[k]
        direction = (indices_arr[1] <= boundary) + (indices_arr[2] <= boundary)

        neighbor = neighbor_ids[boundary]
        side = neighbor_sides[boundary]
        orientation = orientations[boundary]

        u_ll, u_rr = get_surface_node_vars(boundaries_u, equations, boundary)
        u_inner = isequal(side, 1) * u_ll + (1 - isequal(side, 1)) * u_rr
        x = get_node_coords(node_coordinates, equations, boundary)

        boundary_flux_node = boundary_stable_helper(
            boundary_conditions,
            u_inner,
            orientation,
            direction,
            x,
            t,
            surface_flux,
            equations,
        )

        @inbounds begin
            for ii in axes(surface_flux_values, 1)
                surface_flux_values[ii, direction, neighbor] = boundary_flux_node[ii]
            end
        end
    end

    return nothing
end

# Assert 
function cuda_boundary_flux!(
    t,
    mesh::TreeMesh{1},
    boundary_condition::BoundaryConditionPeriodic,
    equations,
    dg::DGSEM,
    cache,
)

    @assert iszero(length(cache.boundaries.orientations))
end

# Launch CUDA kernels to calculate boundary fluxes
function cuda_boundary_flux!(
    t,
    mesh::TreeMesh{1},
    boundary_conditions::NamedTuple,
    equations,
    dg::DGSEM,
    cache,
)

    surface_flux = dg.surface_integral.surface_flux

    n_boundaries_per_direction = CuArray{Int}(cache.boundaries.n_boundaries_per_direction)
    neighbor_ids = CuArray{Int}(cache.boundaries.neighbor_ids)
    neighbor_sides = CuArray{Int}(cache.boundaries.neighbor_sides)
    orientations = CuArray{Int}(cache.boundaries.orientations)
    boundaries_u = CuArray{Float64}(cache.boundaries.u)
    node_coordinates = CuArray{Float64}(cache.boundaries.node_coordinates)
    surface_flux_values = CuArray{Float64}(cache.elements.surface_flux_values)

    lasts = CUDA.zeros(Int, length(n_boundaries_per_direction))
    firsts = CUDA.zeros(Int, length(n_boundaries_per_direction))

    last_first_indices_kernel = @cuda launch = false last_first_indices_kernel!(
        lasts,
        firsts,
        n_boundaries_per_direction,
    )
    last_first_indices_kernel(
        lasts,
        firsts,
        n_boundaries_per_direction;
        configurator_1d(last_first_indices_kernel, lasts)...,
    )

    lasts, firsts = Array(lasts), Array(firsts)
    boundary_arr = CuArray{Int}(firsts[1]:lasts[2])
    indices_arr = CuArray{Int}([firsts[1], firsts[2]])

    size_arr = CuArray{Float64}(undef, length(boundary_arr))

    boundary_flux_kernel = @cuda launch = false boundary_flux_kernel!(
        surface_flux_values,
        boundaries_u,
        node_coordinates,
        t,
        boundary_arr,
        indices_arr,
        neighbor_ids,
        neighbor_sides,
        orientations,
        boundary_conditions,
        equations,
        surface_flux,
    )
    boundary_flux_kernel(
        surface_flux_values,
        boundaries_u,
        node_coordinates,
        t,
        boundary_arr,
        indices_arr,
        neighbor_ids,
        neighbor_sides,
        orientations,
        boundary_conditions,
        equations,
        surface_flux;
        configurator_1d(boundary_flux_kernel, size_arr)...,
    )

    cache.elements.surface_flux_values = surface_flux_values # Automatically copy back to CPU

    return nothing
end

# CUDA kernel for calculating surface integrals along axis x 
function surface_integral_kernel!(du, factor_arr, surface_flux_values)

    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    k = (blockIdx().z - 1) * blockDim().z + threadIdx().z

    if (i <= size(du, 1) && j <= size(du, 2) && k <= size(du, 3))
        @inbounds begin
            du[i, j, k] -= surface_flux_values[i, 1, k] * isequal(j, 1) * factor_arr[1]
            du[i, j, k] +=
                surface_flux_values[i, 2, k] * isequal(j, size(du, 2)) * factor_arr[2]
        end
    end

    return nothing
end

# Launch CUDA kernel to calculate surface integrals
function cuda_surface_integral!(du, mesh::TreeMesh{1}, dg::DGSEM, cache)

    factor_arr = CuArray{Float64}([
        dg.basis.boundary_interpolation[1, 1],
        dg.basis.boundary_interpolation[size(du, 2), 2],
    ])
    surface_flux_values = CuArray{Float64}(cache.elements.surface_flux_values)

    size_arr = CuArray{Float64}(undef, size(du, 1), size(du, 2), size(du, 3))

    surface_integral_kernel =
        @cuda launch = false surface_integral_kernel!(du, factor_arr, surface_flux_values)
    surface_integral_kernel(
        du,
        factor_arr,
        surface_flux_values;
        configurator_3d(surface_integral_kernel, size_arr)...,
    )

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

# Launch CUDA kernel to apply Jacobian to reference element
function cuda_jacobian!(du, mesh::TreeMesh{1}, cache)

    inverse_jacobian = CuArray{Float64}(cache.elements.inverse_jacobian)

    jacobian_kernel = @cuda launch = false jacobian_kernel!(du, inverse_jacobian)
    jacobian_kernel(du, inverse_jacobian; configurator_3d(jacobian_kernel, du)...)

    return nothing
end

# CUDA kernel for calculating source terms
function source_terms_kernel!(
    du,
    u,
    node_coordinates,
    t,
    equations::AbstractEquations{1},
    source_terms::Function,
)

    j = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    k = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if (j <= size(du, 2) && k <= size(du, 3))
        u_local = get_nodes_vars(u, equations, j, k)
        x_local = get_node_coords(node_coordinates, equations, j, k)

        source_terms_node = source_terms(u_local, x_local, t, equations)

        @inbounds begin
            for ii in axes(du, 1)
                du[ii, j, k] += source_terms_node[ii]
            end
        end
    end

    return nothing
end

# Return nothing             
function cuda_sources!(
    du,
    u,
    t,
    source_terms::Nothing,
    equations::AbstractEquations{1},
    cache,
)

    return nothing
end

# Launch CUDA kernel to calculate source terms 
function cuda_sources!(du, u, t, source_terms, equations::AbstractEquations{1}, cache)

    node_coordinates = CuArray{Float64}(cache.elements.node_coordinates)

    size_arr = CuArray{Float64}(undef, size(du, 2), size(du, 3))

    source_terms_kernel = @cuda launch = false source_terms_kernel!(
        du,
        u,
        node_coordinates,
        t,
        equations,
        source_terms,
    )
    source_terms_kernel(
        du,
        u,
        node_coordinates,
        t,
        equations,
        source_terms;
        configurator_2d(source_terms_kernel, size_arr)...,
    )

    return nothing
end

# Pack kernels into `rhs!()`
#################################################################################

# For tests
#################################################################################
#= du, u = copy_to_gpu!(du, u)

cuda_volume_integral!(
    du, u, mesh,
    have_nonconservative_terms(equations), equations,
    solver.volume_integral, solver)

cuda_prolong2interfaces!(u, mesh, cache)

cuda_interface_flux!(
    mesh, have_nonconservative_terms(equations),
    equations, solver, cache)

cuda_prolong2boundaries!(u, mesh,
    boundary_conditions, cache)

cuda_boundary_flux!(t, mesh, boundary_conditions,
    equations, solver, cache)

cuda_surface_integral!(du, mesh, solver, cache)

cuda_jacobian!(du, mesh, cache)

cuda_sources!(du, u, t,
    source_terms, equations, cache)

du, u = copy_to_cpu!(du, u) =#



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

prolong2boundaries!(cache, u, mesh, equations,
    solver.surface_integral, solver)

calc_boundary_flux!(cache, t, boundary_conditions, mesh, equations,
    solver.surface_integral, solver)

calc_surface_integral!(
    du, u, mesh, equations, solver.surface_integral, solver, cache)

apply_jacobian!(du, mesh, equations, solver, cache)

calc_sources!(du, u, t,
    source_terms, equations, solver, cache) =#
