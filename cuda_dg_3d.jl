# Remove it after first run to avoid recompilation
#= include("header.jl") =#

# Set random seed for random tests
Random.seed!(123)

# Use the target test header file
#= include("test/advection_basic_3d.jl") =#
#= include("test/euler_ec_3d.jl") =#
#= include("test/euler_source_terms_3d.jl") =#
include("test/hypdiff_nonperiodic_3d.jl")

# Kernel configurators 
#################################################################################

# CUDA kernel configurator for 1D array computing
function configurator_1d(kernel::CUDA.HostKernel, array::CuArray{<:Any,1})
    config = launch_configuration(kernel.fun)

    threads = min(length(array), config.threads)
    blocks = cld(length(array), threads)

    return (threads=threads, blocks=blocks)
end

# CUDA kernel configurator for 2D array computing
function configurator_2d(kernel::CUDA.HostKernel, array::CuArray{<:Any,2})
    config = launch_configuration(kernel.fun)

    threads = Tuple(fill(Int(floor((min(maximum(size(array)), config.threads))^(1 / 2))), 2))
    blocks = map(cld, size(array), threads)

    return (threads=threads, blocks=blocks)
end

# CUDA kernel configurator for 3D array computing
function configurator_3d(kernel::CUDA.HostKernel, array::CuArray{<:Any,3})
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

# CUDA kernel for calculating fluxes along normal direction 1, 2, 3
function flux_kernel!(flux_arr1, flux_arr2, flux_arr3, u,
    equations::AbstractEquations{3}, flux::Function)

    j = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    k = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if (j <= size(u, 2)^3 && k <= size(u, 5))
        j1 = div(j - 1, size(u, 2)^2) + 1
        j2 = div(rem(j - 1, size(u, 2)^2), size(u, 2)) + 1
        j3 = rem(rem(j - 1, size(u, 2)^2), size(u, 2)) + 1

        u_node = get_nodes_vars(u, equations, j1, j2, j3, k)

        flux_node1 = flux(u_node, 1, equations)
        flux_node2 = flux(u_node, 2, equations)
        flux_node3 = flux(u_node, 3, equations)

        @inbounds begin
            for ii in axes(u, 1)
                flux_arr1[ii, j1, j2, j3, k] = flux_node1[ii]
                flux_arr2[ii, j1, j2, j3, k] = flux_node2[ii]
                flux_arr3[ii, j1, j2, j3, k] = flux_node3[ii]
            end
        end
    end

    return nothing
end

# CUDA kernel for calculating weak form
function weak_form_kernel!(du, derivative_dhat, flux_arr1, flux_arr2, flux_arr3)

    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    k = (blockIdx().z - 1) * blockDim().z + threadIdx().z

    if (i <= size(du, 1) && j <= size(du, 2)^3 && k <= size(du, 5))
        j1 = div(j - 1, size(du, 2)^2) + 1
        j2 = div(rem(j - 1, size(du, 2)^2), size(du, 2)) + 1
        j3 = rem(rem(j - 1, size(du, 2)^2), size(du, 2)) + 1

        @inbounds begin
            for ii in axes(du, 2)
                du[i, j1, j2, j3, k] += derivative_dhat[j1, ii] * flux_arr1[i, ii, j2, j3, k]
                du[i, j1, j2, j3, k] += derivative_dhat[j2, ii] * flux_arr2[i, j1, ii, j3, k]
                du[i, j1, j2, j3, k] += derivative_dhat[j3, ii] * flux_arr3[i, j1, j2, ii, k]
            end
        end
    end

    return nothing
end

# CUDA kernel for calculating volume fluxes in direction x, y, z
function volume_flux_kernel!(volume_flux_arr1, volume_flux_arr2, volume_flux_arr3, u,
    equations::AbstractEquations{3}, volume_flux::Function)

    j = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    k = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if (j <= size(u, 2)^4 && k <= size(u, 5))
        j1 = div(j - 1, size(u, 2)^3) + 1
        j2 = div(rem(j - 1, size(u, 2)^3), size(u, 2)^2) + 1
        j3 = div(rem(j - 1, size(u, 2)^2), size(u, 2)) + 1
        j4 = rem(j - 1, size(u, 2)) + 1

        u_node = get_nodes_vars(u, equations, j1, j2, j3, k)
        u_node1 = get_nodes_vars(u, equations, j4, j2, j3, k)
        u_node2 = get_nodes_vars(u, equations, j1, j4, j3, k)
        u_node3 = get_nodes_vars(u, equations, j1, j2, j4, k)

        volume_flux_node1 = volume_flux(u_node, u_node1, 1, equations)
        volume_flux_node2 = volume_flux(u_node, u_node2, 2, equations)
        volume_flux_node3 = volume_flux(u_node, u_node3, 3, equations)

        @inbounds begin
            for ii in axes(u, 1)
                volume_flux_arr1[ii, j1, j4, j2, j3, k] = volume_flux_node1[ii]
                volume_flux_arr2[ii, j1, j2, j4, j3, k] = volume_flux_node2[ii]
                volume_flux_arr3[ii, j1, j2, j3, j4, k] = volume_flux_node3[ii]
            end
        end
    end

    return nothing
end

# CUDA kernel for calculating volume integrals
function volume_integral_kernel!(du, derivative_split, volume_flux_arr1, volume_flux_arr2, volume_flux_arr3)

    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    k = (blockIdx().z - 1) * blockDim().z + threadIdx().z

    if (i <= size(du, 1) && j <= size(du, 2)^3 && k <= size(du, 5))
        j1 = div(j - 1, size(du, 2)^2) + 1
        j2 = div(rem(j - 1, size(du, 2)^2), size(du, 2)) + 1
        j3 = rem(rem(j - 1, size(du, 2)^2), size(du, 2)) + 1

        @inbounds begin
            for ii in axes(du, 2)
                du[i, j1, j2, j3, k] += derivative_split[j1, ii] * volume_flux_arr1[i, j1, ii, j2, j3, k]
                du[i, j1, j2, j3, k] += derivative_split[j2, ii] * volume_flux_arr2[i, j1, j2, ii, j3, k]
                du[i, j1, j2, j3, k] += derivative_split[j3, ii] * volume_flux_arr3[i, j1, j2, j3, ii, k]
            end
        end
    end

    return nothing
end

# Launch CUDA kernels to calculate volume integrals
function cuda_volume_integral!(du, u, mesh::TreeMesh{3},
    nonconservative_terms, equations,
    volume_integral::VolumeIntegralWeakForm, dg::DGSEM)

    derivative_dhat = CuArray{Float32}(dg.basis.derivative_dhat)
    flux_arr1 = similar(u)
    flux_arr2 = similar(u)
    flux_arr3 = similar(u)

    size_arr = CuArray{Float32}(undef, size(u, 2)^3, size(u, 5))

    flux_kernel = @cuda launch = false flux_kernel!(flux_arr1, flux_arr2, flux_arr3, u, equations, flux)
    flux_kernel(flux_arr1, flux_arr2, flux_arr3, u, equations; configurator_2d(flux_kernel, size_arr)...)

    size_arr = CuArray{Float32}(undef, size(du, 1), size(du, 2)^3, size(du, 5))

    weak_form_kernel = @cuda launch = false weak_form_kernel!(du, derivative_dhat, flux_arr1, flux_arr2, flux_arr3)
    weak_form_kernel(du, derivative_dhat, flux_arr1, flux_arr2, flux_arr3; configurator_3d(weak_form_kernel, size_arr)...)

    return nothing
end

# Launch CUDA kernels to calculate volume integrals
function cuda_volume_integral!(du, u, mesh::TreeMesh{3},
    nonconservative_terms::False, equations,
    volume_integral::VolumeIntegralFluxDifferencing, dg::DGSEM)

    volume_flux = volume_integral.volume_flux

    derivative_split = CuArray{Float32}(dg.basis.derivative_split)
    volume_flux_arr1 = CuArray{Float32}(undef, size(u, 1), size(u, 2), size(u, 2), size(u, 2), size(u, 2), size(u, 5))
    volume_flux_arr2 = CuArray{Float32}(undef, size(u, 1), size(u, 2), size(u, 2), size(u, 2), size(u, 2), size(u, 5))
    volume_flux_arr3 = CuArray{Float32}(undef, size(u, 1), size(u, 2), size(u, 2), size(u, 2), size(u, 2), size(u, 5))

    size_arr = CuArray{Float32}(undef, size(u, 2)^4, size(u, 5))

    volume_flux_kernel = @cuda launch = false volume_flux_kernel!(volume_flux_arr1, volume_flux_arr2, volume_flux_arr3, u, equations, volume_flux)
    volume_flux_kernel(volume_flux_arr1, volume_flux_arr2, volume_flux_arr3, u, equations, volume_flux; configurator_2d(volume_flux_kernel, size_arr)...)

    size_arr = CuArray{Float32}(undef, size(du, 1), size(du, 2)^3, size(du, 5))

    volume_integral_kernel = @cuda launch = false volume_integral_kernel!(du, derivative_split, volume_flux_arr1, volume_flux_arr2, volume_flux_arr3)
    volume_integral_kernel(du, derivative_split, volume_flux_arr1, volume_flux_arr2, volume_flux_arr3; configurator_3d(volume_integral_kernel, size_arr)...)

    return nothing
end

# CUDA kernel for prolonging two interfaces in direction x, y, z
function prolong_interfaces_kernel!(interfaces_u, u, neighbor_ids, orientations)

    j = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    k = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if (j <= size(interfaces_u, 2) * size(interfaces_u, 3)^2 && k <= size(interfaces_u, 5))
        j1 = div(j - 1, size(interfaces_u, 3)^2) + 1
        j2 = div(rem(j - 1, size(interfaces_u, 3)^2), size(interfaces_u, 3)) + 1
        j3 = rem(rem(j - 1, size(interfaces_u, 3)^2), size(interfaces_u, 3)) + 1

        orientation = orientations[k]
        left_element = neighbor_ids[1, k]
        right_element = neighbor_ids[2, k]

        @inbounds begin
            interfaces_u[1, j1, j2, j3, k] = u[j1,
                isequal(orientation, 1)*size(u, 2)+isequal(orientation, 2)*j2+isequal(orientation, 3)*j2,
                isequal(orientation, 1)*j2+isequal(orientation, 2)*size(u, 2)+isequal(orientation, 3)*j3,
                isequal(orientation, 1)*j3+isequal(orientation, 2)*j3+isequal(orientation, 3)*size(u, 2),
                left_element]
            interfaces_u[2, j1, j2, j3, k] = u[j1,
                isequal(orientation, 1)*1+isequal(orientation, 2)*j2+isequal(orientation, 3)*j2,
                isequal(orientation, 1)*j2+isequal(orientation, 2)*1+isequal(orientation, 3)*j3,
                isequal(orientation, 1)*j3+isequal(orientation, 2)*j3+isequal(orientation, 3)*1,
                right_element]
        end
    end

    return nothing
end

# Launch CUDA kernel to prolong solution to interfaces
function cuda_prolong2interfaces!(u, mesh::TreeMesh{3}, cache)

    neighbor_ids = CuArray{Int}(cache.interfaces.neighbor_ids)
    orientations = CuArray{Int}(cache.interfaces.orientations)
    interfaces_u = CuArray{Float32}(cache.interfaces.u)

    size_arr = CuArray{Float32}(undef, size(interfaces_u, 2) * size(interfaces_u, 3)^2, size(interfaces_u, 5))

    prolong_interfaces_kernel = @cuda launch = false prolong_interfaces_kernel!(interfaces_u, u, neighbor_ids, orientations)
    prolong_interfaces_kernel(interfaces_u, u, neighbor_ids, orientations; configurator_2d(prolong_interfaces_kernel, size_arr)...)

    cache.interfaces.u = interfaces_u  # Automatically copy back to CPU

    return nothing
end

# CUDA kernel for calculating surface fluxes 
function surface_flux_kernel!(surface_flux_arr, interfaces_u, orientations,
    equations::AbstractEquations{3}, surface_flux::Any)

    j2 = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j3 = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    k = (blockIdx().z - 1) * blockDim().z + threadIdx().z

    if (j2 <= size(surface_flux_arr, 3) && j3 <= size(surface_flux_arr, 4) && k <= size(surface_flux_arr, 5))
        u_ll, u_rr = get_surface_node_vars(interfaces_u, equations, j2, j3, k)
        orientation = orientations[k]

        surface_flux_node = surface_flux(u_ll, u_rr, orientation, equations)

        @inbounds begin
            for j1j1 in axes(surface_flux_arr, 2)
                surface_flux_arr[1, j1j1, j2, j3, k] = surface_flux_node[j1j1]
            end
        end
    end

    return nothing
end

# CUDA kernel for setting interface fluxes on orientation 1, 2, 3
function interface_flux_kernel!(surface_flux_values, surface_flux_arr, neighbor_ids, orientations)

    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    k = (blockIdx().z - 1) * blockDim().z + threadIdx().z

    if (i <= size(surface_flux_values, 1) && j <= size(surface_flux_arr, 3)^2 && k <= size(surface_flux_arr, 5))
        j1 = div(j - 1, size(surface_flux_arr, 3)) + 1
        j2 = rem(j - 1, size(surface_flux_arr, 3)) + 1

        left_id = neighbor_ids[1, k]
        right_id = neighbor_ids[2, k]

        left_direction = 2 * orientations[k]
        right_direction = 2 * orientations[k] - 1

        @inbounds begin
            surface_flux_values[i, j1, j2, left_direction, left_id] = surface_flux_arr[1, i, j1, j2, k]
            surface_flux_values[i, j1, j2, right_direction, right_id] = surface_flux_arr[1, i, j1, j2, k]
        end
    end

    return nothing
end

# Launch CUDA kernels to calculate interface fluxes
function cuda_interface_flux!(mesh::TreeMesh{3}, nonconservative_terms::False,
    equations, dg::DGSEM, cache)

    surface_flux = dg.surface_integral.surface_flux

    neighbor_ids = CuArray{Int}(cache.interfaces.neighbor_ids)
    orientations = CuArray{Int}(cache.interfaces.orientations)
    interfaces_u = CuArray{Float32}(cache.interfaces.u)
    surface_flux_arr = CuArray{Float32}(undef, 1, size(interfaces_u)[2:end]...)
    surface_flux_values = CuArray{Float32}(cache.elements.surface_flux_values)

    size_arr = CuArray{Float32}(undef, size(interfaces_u, 3), size(interfaces_u, 4), size(interfaces_u, 5))

    surface_flux_kernel = @cuda launch = false surface_flux_kernel!(surface_flux_arr, interfaces_u, orientations, equations, surface_flux)
    surface_flux_kernel(surface_flux_arr, interfaces_u, orientations, equations, surface_flux; configurator_3d(surface_flux_kernel, size_arr)...)

    size_arr = CuArray{Float32}(undef, size(surface_flux_values, 1), size(interfaces_u, 3)^2, size(interfaces_u, 5))

    interface_flux_kernel = @cuda launch = false interface_flux_kernel!(surface_flux_values, surface_flux_arr, neighbor_ids, orientations)
    interface_flux_kernel(surface_flux_values, surface_flux_arr, neighbor_ids, orientations; configurator_3d(interface_flux_kernel, size_arr)...)

    cache.elements.surface_flux_values = surface_flux_values # Automatically copy back to CPU

    return nothing
end

# CUDA kernel for prolonging two boundaries in direction x, y, z
function prolong_boundaries_kernel!(boundaries_u, u, neighbor_ids, neighbor_sides, orientations)

    j = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    k = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if (j <= size(boundaries_u, 2) * size(boundaries_u, 3)^2 && k <= size(boundaries_u, 5))
        j1 = div(j - 1, size(boundaries_u, 3)^2) + 1
        j2 = div(rem(j - 1, size(boundaries_u, 3)^2), size(boundaries_u, 3)) + 1
        j3 = rem(rem(j - 1, size(boundaries_u, 3)^2), size(boundaries_u, 3)) + 1

        element = neighbor_ids[k]
        side = neighbor_sides[k]
        orientation = orientations[k]

        @inbounds begin
            boundaries_u[1, j1, j2, j3, k] = u[j1,
                isequal(orientation, 1)*size(u, 2)+isequal(orientation, 2)*j2+isequal(orientation, 3)*j2,
                isequal(orientation, 1)*j2+isequal(orientation, 2)*size(u, 2)+isequal(orientation, 3)*j3,
                isequal(orientation, 1)*j3+isequal(orientation, 2)*j3+isequal(orientation, 3)*size(u, 2),
                element] * isequal(side, 1) # Set to 0 instead of NaN
            boundaries_u[2, j1, j2, j3, k] = u[j1,
                isequal(orientation, 1)*1+isequal(orientation, 2)*j2+isequal(orientation, 3)*j2,
                isequal(orientation, 1)*j2+isequal(orientation, 2)*1+isequal(orientation, 3)*j3,
                isequal(orientation, 1)*j3+isequal(orientation, 2)*j3+isequal(orientation, 3)*1,
                element] * (1 - isequal(side, 1)) # Set to 0 instead of NaN
        end
    end

    return nothing
end

# Assert 
function cuda_prolong2boundaries!(u, mesh::TreeMesh{3}, cache)

    @assert isequal(length(cache.boundaries.orientations), 1)
end

# Launch CUDA kernel to prolong solution to boundaries
function cuda_prolong2boundaries!(u, mesh::TreeMesh{3}, cache)

    neighbor_ids = CuArray{Int}(cache.boundaries.neighbor_ids)
    neighbor_sides = CuArray{Int}(cache.boundaries.neighbor_sides)
    orientations = CuArray{Int}(cache.boundaries.orientations)
    boundaries_u = CuArray{Float32}(cache.boundaries.u)

    size_arr = CuArray{Float32}(undef, size(boundaries_u, 2) * size(boundaries_u, 3)^2, size(boundaries_u, 5))

    prolong_boundaries_kernel = @cuda launch = false prolong_boundaries_kernel!(boundaries_u, u, neighbor_ids, neighbor_sides, orientations)
    prolong_boundaries_kernel(boundaries_u, u, neighbor_ids, neighbor_sides, orientations; configurator_2d(prolong_boundaries_kernel, size_arr)...)

    cache.boundaries.u = boundaries_u  # Automatically copy back to CPU

    return nothing
end

# CUDA kernel for getting last and first indices
function last_first_indices_kernel!(lasts, firsts, n_boundaries_per_direction)

    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    if (i <= length(n_boundaries_per_direction))
        @inbounds begin
            for ii in 1:i
                lasts[i] += n_boundaries_per_direction[ii]
            end
            firsts[i] = lasts[i] - n_boundaries_per_direction[i] + 1
        end
    end

    return nothing
end

# CUDA kernel for calculating boundary fluxes on direction 1, 2, 3, 4, 5, 6
function boundary_flux_kernel!(surface_flux_values, boundaries_u, node_coordinates, t,
    boundary_arr, indices_arr,
    neighbor_ids, neighbor_sides, orientations,
    boundary_conditions::NamedTuple, equations::AbstractEquations{3}, surface_flux::Any)

    j = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    k = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if (j <= size(surface_flux_values, 2)^2 && k <= length(boundary_arr))
        j1 = div(j - 1, size(surface_flux_values, 2)) + 1
        j2 = rem(j - 1, size(surface_flux_values, 2)) + 1

        boundary = boundary_arr[k]
        direction =
            (indices_arr[1] <= boundary) + (indices_arr[2] <= boundary) +
            (indices_arr[3] <= boundary) + (indices_arr[4] <= boundary) +
            (indices_arr[5] <= boundary) + (indices_arr[6] <= boundary)

        neighbor = neighbor_ids[boundary]
        side = neighbor_sides[boundary]
        orientation = orientations[boundary]

        u_ll, u_rr = get_surface_node_vars(boundaries_u, equations, j1, j2, boundary)
        u_inner = isequal(side, 1) * u_ll + (1 - isequal(side, 1)) * u_rr
        x = get_node_coords(node_coordinates, equations, j1, j2, boundary)

        boundary_condition = boundary_conditions[direction]
        boundary_flux_node = boundary_condition(u_inner, orientation, direction, x, t, surface_flux, equations)

        @inbounds begin
            for ii in axes(surface_flux_values, 1)
                surface_flux_values[ii, j1, j2, direction, neighbor] = boundary_flux_node[ii]
            end
        end
    end

    return nothing
end

# Assert 
function cuda_boundary_flux!(t, mesh::TreeMesh{3}, boundary_condition::BoundaryConditionPeriodic,
    equations, dg::DGSEM, cache)

    @assert isequal(length(cache.boundaries.orientations), 1)
end

# Launch CUDA kernels to calculate boundary fluxes
function cuda_boundary_flux!(t, mesh::TreeMesh{3}, boundary_conditions::NamedTuple,
    equations, dg::DGSEM, cache)

    surface_flux = dg.surface_integral.surface_flux

    n_boundaries_per_direction = CuArray{Int}(cache.boundaries.n_boundaries_per_direction)
    neighbor_ids = CuArray{Int}(cache.boundaries.neighbor_ids)
    neighbor_sides = CuArray{Int}(cache.boundaries.neighbor_sides)
    orientations = CuArray{Int}(cache.boundaries.orientations)
    boundaries_u = CuArray{Float32}(cache.boundaries.u)
    node_coordinates = CuArray{Float32}(cache.boundaries.node_coordinates)
    surface_flux_values = CuArray{Float32}(cache.elements.surface_flux_values)

    lasts = CUDA.zeros(Int, length(n_boundaries_per_direction))
    firsts = CUDA.zeros(Int, length(n_boundaries_per_direction))

    last_first_indices_kernel = @cuda launch = false last_first_indices_kernel!(lasts, firsts, n_boundaries_per_direction)
    last_first_indices_kernel(lasts, firsts, n_boundaries_per_direction; configurator_1d(last_first_indices_kernel, lasts)...)

    lasts, firsts = Array(lasts), Array(firsts)
    boundary_arr = CuArray{Int}(firsts[1]:lasts[6])
    indices_arr = CuArray{Int}([firsts[1], firsts[2], firsts[3], firsts[4], firsts[5], firsts[6]])

    size_arr = CuArray{Float32}(undef, size(surface_flux_values, 2)^2, length(boundary_arr))

    boundary_flux_kernel = @cuda launch = false boundary_flux_kernel!(surface_flux_values, boundaries_u, node_coordinates, t,
        boundary_arr, indices_arr, neighbor_ids, neighbor_sides, orientations, boundary_conditions, equations, surface_flux)
    boundary_flux_kernel(surface_flux_values, boundaries_u, node_coordinates, t, boundary_arr, indices_arr, neighbor_ids, neighbor_sides,
        orientations, boundary_conditions, equations, surface_flux; configurator_2d(boundary_flux_kernel, size_arr)...)

    cache.elements.surface_flux_values = surface_flux_values # Automatically copy back to CPU

    return nothing
end

# CUDA kernel for calculating surface integrals along axis x, y, z
function surface_integral_kernel!(du, factor_arr, surface_flux_values)

    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    k = (blockIdx().z - 1) * blockDim().z + threadIdx().z

    if (i <= size(du, 1) && j <= size(du, 2)^3 && k <= size(du, 5))
        j1 = div(j - 1, size(du, 2)^2) + 1
        j2 = div(rem(j - 1, size(du, 2)^2), size(du, 2)) + 1
        j3 = rem(rem(j - 1, size(du, 2)^2), size(du, 2)) + 1

        @inbounds begin
            du[i, j1, j2, j3, k] -= (surface_flux_values[i, j2, j3, 1, k] * isequal(j1, 1)
                                     +
                                     surface_flux_values[i, j1, j3, 3, k] * isequal(j2, 1)
                                     +
                                     surface_flux_values[i, j1, j2, 5, k] * isequal(j3, 1)) * factor_arr[1]
            du[i, j1, j2, j3, k] += (surface_flux_values[i, j2, j3, 2, k] * isequal(j1, size(du, 2))
                                     +
                                     surface_flux_values[i, j1, j3, 4, k] * isequal(j2, size(du, 2))
                                     +
                                     surface_flux_values[i, j1, j2, 6, k] * isequal(j3, size(du, 2))) * factor_arr[2]
        end
    end

    return nothing
end

# Launch CUDA kernel to calculate surface integrals
function cuda_surface_integral!(du, mesh::TreeMesh{3}, dg::DGSEM, cache) # surface_integral

    factor_arr = CuArray{Float32}([dg.basis.boundary_interpolation[1, 1], dg.basis.boundary_interpolation[size(du, 2), 2]])
    surface_flux_values = CuArray{Float32}(cache.elements.surface_flux_values)

    size_arr = CuArray{Float32}(undef, size(du, 1), size(du, 2)^3, size(du, 5))

    surface_integral_kernel = @cuda launch = false surface_integral_kernel!(du, factor_arr, surface_flux_values)
    surface_integral_kernel(du, factor_arr, surface_flux_values; configurator_3d(surface_integral_kernel, size_arr)...)

    return nothing
end

# CUDA kernel for applying inverse Jacobian 
function jacobian_kernel!(du, inverse_jacobian)

    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    k = (blockIdx().z - 1) * blockDim().z + threadIdx().z

    if (i <= size(du, 1) && j <= size(du, 2)^3 && k <= size(du, 5))
        j1 = div(j - 1, size(du, 2)^2) + 1
        j2 = div(rem(j - 1, size(du, 2)^2), size(du, 2)) + 1
        j3 = rem(rem(j - 1, size(du, 2)^2), size(du, 2)) + 1

        @inbounds du[i, j1, j2, j3, k] *= -inverse_jacobian[k]
    end

    return nothing
end

# Launch CUDA kernel to apply Jacobian to reference element
function cuda_jacobian!(du, mesh::TreeMesh{3}, cache)

    inverse_jacobian = CuArray{Float32}(cache.elements.inverse_jacobian)

    size_arr = CuArray{Float32}(undef, size(du, 1), size(du, 2)^3, size(du, 5))

    jacobian_kernel = @cuda launch = false jacobian_kernel!(du, inverse_jacobian)
    jacobian_kernel(du, inverse_jacobian; configurator_3d(jacobian_kernel, size_arr)...)

    return nothing
end

# CUDA kernel for calculating source terms
function source_terms_kernel!(du, u, node_coordinates, t,
    equations::AbstractEquations{3}, source_terms::Function)

    j = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    k = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if (j <= size(du, 2)^3 && k <= size(du, 5))
        j1 = div(j - 1, size(du, 2)^2) + 1
        j2 = div(rem(j - 1, size(du, 2)^2), size(du, 2)) + 1
        j3 = rem(rem(j - 1, size(du, 2)^2), size(du, 2)) + 1

        u_local = get_nodes_vars(u, equations, j1, j2, j3, k)
        x_local = get_node_coords(node_coordinates, equations, j1, j2, j3, k)

        source_terms_node = source_terms(u_local, x_local, t, equations)

        @inbounds begin
            for ii in axes(du, 1)
                du[ii, j1, j2, j3, k] += source_terms_node[ii]
            end
        end
    end

    return nothing
end

# Return nothing           
function cuda_sources!(du, u, t, source_terms::Nothing,
    equations::AbstractEquations{3}, cache)

    return nothing
end

# Launch CUDA kernel to calculate source terms 
function cuda_sources!(du, u, t, source_terms,
    equations::AbstractEquations{3}, cache)

    node_coordinates = CuArray{Float32}(cache.elements.node_coordinates)

    size_arr = CuArray{Float32}(undef, size(u, 2)^3, size(u, 5))

    source_terms_kernel = @cuda launch = false source_terms_kernel!(du, u, node_coordinates, t, equations, source_terms)
    source_terms_kernel(du, u, node_coordinates, t, equations, source_terms; configurator_2d(source_terms_kernel, size_arr)...)

    return nothing
end

# Pack kernels into `rhs!()`
#################################################################################

# For tests
#################################################################################
du, u = copy_to_gpu!(du, u)

cuda_volume_integral!(
    du, u, mesh,
    have_nonconservative_terms(equations), equations,
    solver.volume_integral, solver)

cuda_prolong2interfaces!(u, mesh, cache)

cuda_interface_flux!(
    mesh, have_nonconservative_terms(equations),
    equations, solver, cache,)

cuda_prolong2boundaries!(u, mesh, cache)

cuda_boundary_flux!(t, mesh, boundary_conditions,
    equations, solver, cache)

#= try
    cuda_boundary_flux!(t, mesh, boundary_conditions,
        equations, solver, cache)
catch err
    code_typed(err; interactive=true)
end =#

#= cuda_surface_integral!(du, mesh, solver, cache)

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

calc_surface_integral!(
    du, u, mesh, equations, solver.surface_integral, solver, cache)

apply_jacobian!(du, mesh, equations, solver, cache)

calc_sources!(du, u, t,
    source_terms, equations, solver, cache) =#
