# Remove it after first run to avoid recompilation
#= include("header.jl") =#

# Use the target test header file
#= include("tests/advection_basic_2d.jl") =#
#= include("tests/euler_ec_2d.jl") =#
#= include("tests/euler_source_terms_2d.jl") =#
#= include("tests/hypdiff_nonperiodic_2d.jl") =#
#= include("tests/advection_mortar_2d.jl") =#
#= include("tests/shallowwater_well_balanced_2d.jl") =#

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

# CUDA kernel for calculating fluxes along normal direction 1, 2
function flux_kernel!(
    flux_arr1,
    flux_arr2,
    u,
    equations::AbstractEquations{2},
    flux::Function,
)

    j = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    k = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if (j <= size(u, 2)^2 && k <= size(u, 4))
        j1 = div(j - 1, size(u, 2)) + 1
        j2 = rem(j - 1, size(u, 2)) + 1

        u_node = get_nodes_vars(u, equations, j1, j2, k)

        flux_node1 = flux(u_node, 1, equations)
        flux_node2 = flux(u_node, 2, equations)

        @inbounds begin
            for ii in axes(u, 1)
                flux_arr1[ii, j1, j2, k] = flux_node1[ii]
                flux_arr2[ii, j1, j2, k] = flux_node2[ii]
            end
        end
    end

    return nothing
end

# CUDA kernel for calculating weak form
function weak_form_kernel!(du, derivative_dhat, flux_arr1, flux_arr2)

    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    k = (blockIdx().z - 1) * blockDim().z + threadIdx().z

    if (i <= size(du, 1) && j <= size(du, 2)^2 && k <= size(du, 4))
        j1 = div(j - 1, size(du, 2)) + 1
        j2 = rem(j - 1, size(du, 2)) + 1

        @inbounds begin
            for ii in axes(du, 2)
                du[i, j1, j2, k] += derivative_dhat[j1, ii] * flux_arr1[i, ii, j2, k]
                du[i, j1, j2, k] += derivative_dhat[j2, ii] * flux_arr2[i, j1, ii, k]
            end
        end
    end

    return nothing
end

# CUDA kernel for calculating volume fluxes in direction x, y
function volume_flux_kernel!(
    volume_flux_arr1,
    volume_flux_arr2,
    u,
    equations::AbstractEquations{2},
    volume_flux::Function,
)

    j = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    k = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if (j <= size(u, 2)^3 && k <= size(u, 4))
        j1 = div(j - 1, size(u, 2)^2) + 1
        j2 = div(rem(j - 1, size(u, 2)^2), size(u, 2)) + 1
        j3 = rem(rem(j - 1, size(u, 2)^2), size(u, 2)) + 1

        u_node = get_nodes_vars(u, equations, j1, j2, k)
        u_node1 = get_nodes_vars(u, equations, j3, j2, k)
        u_node2 = get_nodes_vars(u, equations, j1, j3, k)

        volume_flux_node1 = volume_flux(u_node, u_node1, 1, equations)
        volume_flux_node2 = volume_flux(u_node, u_node2, 2, equations)

        @inbounds begin
            for ii in axes(u, 1)
                volume_flux_arr1[ii, j1, j3, j2, k] = volume_flux_node1[ii]
                volume_flux_arr2[ii, j1, j2, j3, k] = volume_flux_node2[ii]
            end
        end
    end

    return nothing
end

# CUDA kernel for calculating symmetric and nonsymmetric fluxes in direction x, y
function symmetric_noncons_flux_kernel!(
    symmetric_flux_arr1,
    symmetric_flux_arr2,
    noncons_flux_arr1,
    noncons_flux_arr2,
    u,
    derivative_split,
    equations::AbstractEquations{2},
    symmetric_flux::Function,
    nonconservative_flux::Function,
)

    j = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    k = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if (j <= size(u, 2)^3 && k <= size(u, 4))
        j1 = div(j - 1, size(u, 2)^2) + 1
        j2 = div(rem(j - 1, size(u, 2)^2), size(u, 2)) + 1
        j3 = rem(rem(j - 1, size(u, 2)^2), size(u, 2)) + 1

        u_node = get_nodes_vars(u, equations, j1, j2, k)
        u_node1 = get_nodes_vars(u, equations, j3, j2, k)
        u_node2 = get_nodes_vars(u, equations, j1, j3, k)

        symmetric_flux_node1 = symmetric_flux(u_node, u_node1, 1, equations)
        symmetric_flux_node2 = symmetric_flux(u_node, u_node2, 2, equations)

        noncons_flux_node1 = nonconservative_flux(u_node, u_node1, 1, equations)
        noncons_flux_node2 = nonconservative_flux(u_node, u_node2, 2, equations)

        @inbounds begin
            for ii in axes(u, 1)
                symmetric_flux_arr1[ii, j1, j3, j2, k] = symmetric_flux_node1[ii]
                symmetric_flux_arr2[ii, j1, j2, j3, k] = symmetric_flux_node2[ii]
                noncons_flux_arr1[ii, j1, j3, j2, k] =
                    noncons_flux_node1[ii] * derivative_split[j1, j3]
                noncons_flux_arr2[ii, j1, j2, j3, k] =
                    noncons_flux_node2[ii] * derivative_split[j2, j3]
            end
        end
    end

    return nothing
end

# CUDA kernel for calculating volume integrals
function volume_integral_kernel!(du, derivative_split, volume_flux_arr1, volume_flux_arr2)

    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    k = (blockIdx().z - 1) * blockDim().z + threadIdx().z

    if (i <= size(du, 1) && j <= size(du, 2)^2 && k <= size(du, 4))
        j1 = div(j - 1, size(du, 2)) + 1
        j2 = rem(j - 1, size(du, 2)) + 1

        @inbounds begin
            for ii in axes(du, 2)
                du[i, j1, j2, k] +=
                    derivative_split[j1, ii] * volume_flux_arr1[i, j1, ii, j2, k]
                du[i, j1, j2, k] +=
                    derivative_split[j2, ii] * volume_flux_arr2[i, j1, j2, ii, k]
            end
        end
    end

    return nothing
end

# CUDA kernel for calculating symmetric and nonsymmetric volume integrals
function symmetric_noncons_integral_kernel!(
    du,
    derivative_split,
    symmetric_flux_arr1,
    symmetric_flux_arr2,
    noncons_flux_arr1,
    noncons_flux_arr2,
)

    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    k = (blockIdx().z - 1) * blockDim().z + threadIdx().z

    if (i <= size(du, 1) && j <= size(du, 2)^2 && k <= size(du, 4))
        j1 = div(j - 1, size(du, 2)) + 1
        j2 = rem(j - 1, size(du, 2)) + 1

        @inbounds begin
            integral_contribution = 0.0

            for ii in axes(du, 2)
                du[i, j1, j2, k] +=
                    derivative_split[j1, ii] * symmetric_flux_arr1[i, j1, ii, j2, k]
                du[i, j1, j2, k] +=
                    derivative_split[j2, ii] * symmetric_flux_arr2[i, j1, j2, ii, k]
                integral_contribution += noncons_flux_arr1[i, j1, ii, j2, k]
                integral_contribution += noncons_flux_arr2[i, j1, j2, ii, k]
            end

            du[i, j1, j2, k] += 0.5 * integral_contribution
        end
    end

    return nothing
end

# Launch CUDA kernels to calculate volume integrals
function cuda_volume_integral!(
    du,
    u,
    mesh::TreeMesh{2},
    nonconservative_terms,
    equations,
    volume_integral::VolumeIntegralWeakForm,
    dg::DGSEM,
)

    derivative_dhat = CuArray{Float64}(dg.basis.derivative_dhat)
    flux_arr1 = similar(u)
    flux_arr2 = similar(u)

    size_arr = CuArray{Float64}(undef, size(u, 2)^2, size(u, 4))

    flux_kernel =
        @cuda launch = false flux_kernel!(flux_arr1, flux_arr2, u, equations, flux)
    flux_kernel(
        flux_arr1,
        flux_arr2,
        u,
        equations;
        configurator_2d(flux_kernel, size_arr)...,
    )

    size_arr = CuArray{Float64}(undef, size(du, 1), size(du, 2)^2, size(du, 4))

    weak_form_kernel =
        @cuda launch = false weak_form_kernel!(du, derivative_dhat, flux_arr1, flux_arr2)
    weak_form_kernel(
        du,
        derivative_dhat,
        flux_arr1,
        flux_arr2;
        configurator_3d(weak_form_kernel, size_arr)...,
    )

    return nothing
end

# Launch CUDA kernels to calculate volume integrals
function cuda_volume_integral!(
    du,
    u,
    mesh::TreeMesh{2},
    nonconservative_terms::False,
    equations,
    volume_integral::VolumeIntegralFluxDifferencing,
    dg::DGSEM,
)

    volume_flux = volume_integral.volume_flux

    derivative_split = CuArray{Float64}(dg.basis.derivative_split)
    volume_flux_arr1 =
        CuArray{Float64}(undef, size(u, 1), size(u, 2), size(u, 2), size(u, 2), size(u, 4))
    volume_flux_arr2 =
        CuArray{Float64}(undef, size(u, 1), size(u, 2), size(u, 2), size(u, 2), size(u, 4))

    size_arr = CuArray{Float64}(undef, size(u, 2)^3, size(u, 4))

    volume_flux_kernel = @cuda launch = false volume_flux_kernel!(
        volume_flux_arr1,
        volume_flux_arr2,
        u,
        equations,
        volume_flux,
    )
    volume_flux_kernel(
        volume_flux_arr1,
        volume_flux_arr2,
        u,
        equations,
        volume_flux;
        configurator_2d(volume_flux_kernel, size_arr)...,
    )

    size_arr = CuArray{Float64}(undef, size(du, 1), size(du, 2)^2, size(du, 4))

    volume_integral_kernel = @cuda launch = false volume_integral_kernel!(
        du,
        derivative_split,
        volume_flux_arr1,
        volume_flux_arr2,
    )
    volume_integral_kernel(
        du,
        derivative_split,
        volume_flux_arr1,
        volume_flux_arr2;
        configurator_3d(volume_integral_kernel, size_arr)...,
    )

    return nothing
end

# Launch CUDA kernels to calculate volume integrals
function cuda_volume_integral!(
    du,
    u,
    mesh::TreeMesh{2},
    nonconservative_terms::True,
    equations,
    volume_integral::VolumeIntegralFluxDifferencing,
    dg::DGSEM,
)

    symmetric_flux, nonconservative_flux = dg.volume_integral.volume_flux

    derivative_split = CuArray{Float64}(dg.basis.derivative_split)
    symmetric_flux_arr1 =
        CuArray{Float64}(undef, size(u, 1), size(u, 2), size(u, 2), size(u, 2), size(u, 4))
    symmetric_flux_arr2 =
        CuArray{Float64}(undef, size(u, 1), size(u, 2), size(u, 2), size(u, 2), size(u, 4))
    noncons_flux_arr1 =
        CuArray{Float64}(undef, size(u, 1), size(u, 2), size(u, 2), size(u, 2), size(u, 4))
    noncons_flux_arr2 =
        CuArray{Float64}(undef, size(u, 1), size(u, 2), size(u, 2), size(u, 2), size(u, 4))

    size_arr = CuArray{Float64}(undef, size(u, 2)^3, size(u, 4))

    symmetric_noncons_flux_kernel = @cuda launch = false symmetric_noncons_flux_kernel!(
        symmetric_flux_arr1,
        symmetric_flux_arr2,
        noncons_flux_arr1,
        noncons_flux_arr2,
        u,
        derivative_split,
        equations,
        symmetric_flux,
        nonconservative_flux,
    )
    symmetric_noncons_flux_kernel(
        symmetric_flux_arr1,
        symmetric_flux_arr2,
        noncons_flux_arr1,
        noncons_flux_arr2,
        u,
        derivative_split,
        equations,
        symmetric_flux,
        nonconservative_flux;
        configurator_2d(symmetric_noncons_flux_kernel, size_arr)...,
    )

    size_arr = CuArray{Float64}(undef, size(du, 1), size(du, 2)^2, size(du, 4))

    symmetric_noncons_integral_kernel =
        @cuda launch = false symmetric_noncons_integral_kernel!(
            du,
            derivative_split,
            symmetric_flux_arr1,
            symmetric_flux_arr2,
            noncons_flux_arr1,
            noncons_flux_arr2,
        )
    symmetric_noncons_integral_kernel(
        du,
        derivative_split,
        symmetric_flux_arr1,
        symmetric_flux_arr2,
        noncons_flux_arr1,
        noncons_flux_arr2;
        configurator_3d(symmetric_noncons_integral_kernel, size_arr)...,
    )

    return nothing
end

# CUDA kernel for prolonging two interfaces in direction x, y
function prolong_interfaces_kernel!(interfaces_u, u, neighbor_ids, orientations)

    j = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    k = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if (j <= size(interfaces_u, 2) * size(interfaces_u, 3) && k <= size(interfaces_u, 4))
        j1 = div(j - 1, size(interfaces_u, 3)) + 1
        j2 = rem(j - 1, size(interfaces_u, 3)) + 1

        orientation = orientations[k]
        left_element = neighbor_ids[1, k]
        right_element = neighbor_ids[2, k]

        @inbounds begin
            interfaces_u[1, j1, j2, k] = u[
                j1,
                isequal(orientation, 1)*size(u, 2)+isequal(orientation, 2)*j2,
                isequal(orientation, 1)*j2+isequal(orientation, 2)*size(u, 2),
                left_element,
            ]
            interfaces_u[2, j1, j2, k] = u[
                j1,
                isequal(orientation, 1)*1+isequal(orientation, 2)*j2,
                isequal(orientation, 1)*j2+isequal(orientation, 2)*1,
                right_element,
            ]
        end
    end

    return nothing
end

# Launch CUDA kernel to prolong solution to interfaces
function cuda_prolong2interfaces!(u, mesh::TreeMesh{2}, cache)

    neighbor_ids = CuArray{Int}(cache.interfaces.neighbor_ids)
    orientations = CuArray{Int}(cache.interfaces.orientations)
    interfaces_u = CuArray{Float64}(cache.interfaces.u)

    size_arr = CuArray{Float64}(
        undef,
        size(interfaces_u, 2) * size(interfaces_u, 3),
        size(interfaces_u, 4),
    )

    prolong_interfaces_kernel = @cuda launch = false prolong_interfaces_kernel!(
        interfaces_u,
        u,
        neighbor_ids,
        orientations,
    )
    prolong_interfaces_kernel(
        interfaces_u,
        u,
        neighbor_ids,
        orientations;
        configurator_2d(prolong_interfaces_kernel, size_arr)...,
    )

    cache.interfaces.u = interfaces_u  # Automatically copy back to CPU

    return nothing
end

# CUDA kernel for calculating surface fluxes 
function surface_flux_kernel!(
    surface_flux_arr,
    interfaces_u,
    orientations,
    equations::AbstractEquations{2},
    surface_flux::Any,
)

    j2 = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    k = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if (j2 <= size(surface_flux_arr, 3) && k <= size(surface_flux_arr, 4))
        u_ll, u_rr = get_surface_node_vars(interfaces_u, equations, j2, k)
        orientation = orientations[k]

        surface_flux_node = surface_flux(u_ll, u_rr, orientation, equations)

        @inbounds begin
            for j1j1 in axes(surface_flux_arr, 2)
                surface_flux_arr[1, j1j1, j2, k] = surface_flux_node[j1j1]
            end
        end
    end

    return nothing
end

# CUDA kernel for setting interface fluxes on orientation 1, 2
function interface_flux_kernel!(
    surface_flux_values,
    surface_flux_arr,
    neighbor_ids,
    orientations,
)

    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    k = (blockIdx().z - 1) * blockDim().z + threadIdx().z

    if (
        i <= size(surface_flux_values, 1) &&
        j <= size(surface_flux_arr, 3) &&
        k <= size(surface_flux_arr, 4)
    )
        left_id = neighbor_ids[1, k]
        right_id = neighbor_ids[2, k]

        left_direction = 2 * orientations[k]
        right_direction = 2 * orientations[k] - 1

        @inbounds begin
            surface_flux_values[i, j, left_direction, left_id] =
                surface_flux_arr[1, i, j, k]
            surface_flux_values[i, j, right_direction, right_id] =
                surface_flux_arr[1, i, j, k]
        end
    end

    return nothing
end

# Launch CUDA kernels to calculate interface fluxes
function cuda_interface_flux!(
    mesh::TreeMesh{2},
    nonconservative_terms::False,
    equations,
    dg::DGSEM,
    cache,
)

    surface_flux = dg.surface_integral.surface_flux

    neighbor_ids = CuArray{Int}(cache.interfaces.neighbor_ids)
    orientations = CuArray{Int}(cache.interfaces.orientations)
    interfaces_u = CuArray{Float64}(cache.interfaces.u)
    surface_flux_arr = CuArray{Float64}(undef, 1, size(interfaces_u)[2:end]...)
    surface_flux_values = CuArray{Float64}(cache.elements.surface_flux_values)

    size_arr = CuArray{Float64}(undef, size(interfaces_u, 3), size(interfaces_u, 4))

    surface_flux_kernel = @cuda launch = false surface_flux_kernel!(
        surface_flux_arr,
        interfaces_u,
        orientations,
        equations,
        surface_flux,
    )
    surface_flux_kernel(
        surface_flux_arr,
        interfaces_u,
        orientations,
        equations,
        surface_flux;
        configurator_2d(surface_flux_kernel, size_arr)...,
    )

    size_arr = CuArray{Float64}(
        undef,
        size(surface_flux_values, 1),
        size(interfaces_u, 3),
        size(interfaces_u, 4),
    )

    interface_flux_kernel = @cuda launch = false interface_flux_kernel!(
        surface_flux_values,
        surface_flux_arr,
        neighbor_ids,
        orientations,
    )
    interface_flux_kernel(
        surface_flux_values,
        surface_flux_arr,
        neighbor_ids,
        orientations;
        configurator_3d(interface_flux_kernel, size_arr)...,
    )

    cache.elements.surface_flux_values = surface_flux_values # Automatically copy back to CPU

    return nothing
end

# CUDA kernel for prolonging two boundaries in direction x, y
function prolong_boundaries_kernel!(
    boundaries_u,
    u,
    neighbor_ids,
    neighbor_sides,
    orientations,
)

    j = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    k = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if (j <= size(boundaries_u, 2) * size(boundaries_u, 3) && k <= size(boundaries_u, 4))
        j1 = div(j - 1, size(boundaries_u, 3)) + 1
        j2 = rem(j - 1, size(boundaries_u, 3)) + 1

        element = neighbor_ids[k]
        side = neighbor_sides[k]
        orientation = orientations[k]

        @inbounds begin
            boundaries_u[1, j1, j2, k] =
                u[
                    j1,
                    isequal(orientation, 1)*size(u, 2)+isequal(orientation, 2)*j2,
                    isequal(orientation, 1)*j2+isequal(orientation, 2)*size(u, 2),
                    element,
                ] * isequal(side, 1) # Set to 0 instead of NaN
            boundaries_u[2, j1, j2, k] =
                u[
                    j1,
                    isequal(orientation, 1)*1+isequal(orientation, 2)*j2,
                    isequal(orientation, 1)*j2+isequal(orientation, 2)*1,
                    element,
                ] * (1 - isequal(side, 1)) # Set to 0 instead of NaN
        end
    end

    return nothing
end

# Assert 
function cuda_prolong2boundaries!(
    u,
    mesh::TreeMesh{2},
    boundary_condition::BoundaryConditionPeriodic,
    cache,
)

    @assert iszero(length(cache.boundaries.orientations))
end

# Launch CUDA kernel to prolong solution to boundaries
function cuda_prolong2boundaries!(
    u,
    mesh::TreeMesh{2},
    boundary_conditions::NamedTuple,
    cache,
)

    neighbor_ids = CuArray{Int}(cache.boundaries.neighbor_ids)
    neighbor_sides = CuArray{Int}(cache.boundaries.neighbor_sides)
    orientations = CuArray{Int}(cache.boundaries.orientations)
    boundaries_u = CuArray{Float64}(cache.boundaries.u)

    size_arr = CuArray{Float64}(
        undef,
        size(boundaries_u, 2) * size(boundaries_u, 3),
        size(boundaries_u, 4),
    )

    prolong_boundaries_kernel = @cuda launch = false prolong_boundaries_kernel!(
        boundaries_u,
        u,
        neighbor_ids,
        neighbor_sides,
        orientations,
    )
    prolong_boundaries_kernel(
        boundaries_u,
        u,
        neighbor_ids,
        neighbor_sides,
        orientations;
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

# CUDA kernel for calculating boundary fluxes in direction 1, 2, 3, 4
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
    equations::AbstractEquations{2},
    surface_flux::Any,
)

    j = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    k = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if (j <= size(surface_flux_values, 2) && k <= length(boundary_arr))
        boundary = boundary_arr[k]
        direction =
            (indices_arr[1] <= boundary) +
            (indices_arr[2] <= boundary) +
            (indices_arr[3] <= boundary) +
            (indices_arr[4] <= boundary)

        neighbor = neighbor_ids[boundary]
        side = neighbor_sides[boundary]
        orientation = orientations[boundary]

        u_ll, u_rr = get_surface_node_vars(boundaries_u, equations, j, boundary)
        u_inner = isequal(side, 1) * u_ll + (1 - isequal(side, 1)) * u_rr
        x = get_node_coords(node_coordinates, equations, j, boundary)

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
                surface_flux_values[ii, j, direction, neighbor] = boundary_flux_node[ii]
            end
        end
    end

    return nothing
end

# Assert 
function cuda_boundary_flux!(
    t,
    mesh::TreeMesh{2},
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
    mesh::TreeMesh{2},
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
    boundary_arr = CuArray{Int}(firsts[1]:lasts[4])
    indices_arr = CuArray{Int}([firsts[1], firsts[2], firsts[3], firsts[4]])

    size_arr = CuArray{Float64}(undef, size(surface_flux_values, 2), length(boundary_arr))

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
        configurator_2d(boundary_flux_kernel, size_arr)...,
    )

    cache.elements.surface_flux_values = surface_flux_values # Automatically copy back to CPU

    return nothing
end

# CUDA kernel for copying data small to small on mortars
function prolong_mortars_small2small_kernel!(
    u_upper,
    u_lower,
    u,
    neighbor_ids,
    large_sides,
    orientations,
)

    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    k = (blockIdx().z - 1) * blockDim().z + threadIdx().z

    if (i <= size(u_upper, 2) && j <= size(u_upper, 3) && k <= size(u_upper, 4))
        large_side = large_sides[k]
        orientation = orientations[k]

        lower_element = neighbor_ids[1, k]
        upper_element = neighbor_ids[2, k]

        @inbounds begin
            u_upper[2, i, j, k] =
                u[
                    i,
                    isequal(orientation, 1)*1+isequal(orientation, 2)*j,
                    isequal(orientation, 1)*j+isequal(orientation, 2)*1,
                    upper_element,
                ] * isequal(large_side, 1)

            u_lower[2, i, j, k] =
                u[
                    i,
                    isequal(orientation, 1)*1+isequal(orientation, 2)*j,
                    isequal(orientation, 1)*j+isequal(orientation, 2)*1,
                    lower_element,
                ] * isequal(large_side, 1)

            u_upper[1, i, j, k] =
                u[
                    i,
                    isequal(orientation, 1)*size(u, 2)+isequal(orientation, 2)*j,
                    isequal(orientation, 1)*j+isequal(orientation, 2)*size(u, 2),
                    upper_element,
                ] * isequal(large_side, 2)

            u_lower[1, i, j, k] =
                u[
                    i,
                    isequal(orientation, 1)*size(u, 2)+isequal(orientation, 2)*j,
                    isequal(orientation, 1)*j+isequal(orientation, 2)*size(u, 2),
                    lower_element,
                ] * isequal(large_side, 2)
        end
    end

    return nothing
end

# CUDA kernel for interpolating data large to small on mortars
function prolong_mortars_large2small_kernel!(
    u_upper,
    u_lower,
    u,
    forward_upper,
    forward_lower,
    neighbor_ids,
    large_sides,
    orientations,
)

    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    k = (blockIdx().z - 1) * blockDim().z + threadIdx().z

    if (i <= size(u_upper, 2) && j <= size(u_upper, 3) && k <= size(u_upper, 4))
        large_side = large_sides[k]
        orientation = orientations[k]
        large_element = neighbor_ids[3, k]

        leftright = large_side

        @inbounds begin
            for ii in axes(forward_upper, 2)
                u_upper[leftright, i, j, k] +=
                    forward_upper[j, ii] *
                    u[
                        i,
                        isequal(orientation, 1)*size(u, 2)+isequal(orientation, 2)*ii,
                        isequal(orientation, 1)*ii+isequal(orientation, 2)*size(u, 2),
                        large_element,
                    ] *
                    isequal(large_side, 1)
                u_lower[leftright, i, j, k] +=
                    forward_lower[j, ii] *
                    u[
                        i,
                        isequal(orientation, 1)*size(u, 2)+isequal(orientation, 2)*ii,
                        isequal(orientation, 1)*ii+isequal(orientation, 2)*size(u, 2),
                        large_element,
                    ] *
                    isequal(large_side, 1)
            end

            for ii in axes(forward_lower, 2)
                u_upper[leftright, i, j, k] +=
                    forward_upper[j, ii] *
                    u[
                        i,
                        isequal(orientation, 1)*1+isequal(orientation, 2)*ii,
                        isequal(orientation, 1)*ii+isequal(orientation, 2)*1,
                        large_element,
                    ] *
                    isequal(large_side, 2)
                u_lower[leftright, i, j, k] +=
                    forward_lower[j, ii] *
                    u[
                        i,
                        isequal(orientation, 1)*1+isequal(orientation, 2)*ii,
                        isequal(orientation, 1)*ii+isequal(orientation, 2)*1,
                        large_element,
                    ] *
                    isequal(large_side, 2)
            end
        end
    end

    return nothing
end

# Assert 
function cuda_prolong2mortars!(u, mesh::TreeMesh{2}, dg::DGSEM, cache_mortars::True, cache)

    @assert iszero(length(cache.mortars.orientations))
end

# Launch CUDA kernels to prolong solution to mortars
function cuda_prolong2mortars!(u, mesh::TreeMesh{2}, dg::DGSEM, cache_mortars::False, cache)

    neighbor_ids = CuArray{Int}(cache.mortars.neighbor_ids)
    large_sides = CuArray{Int}(cache.mortars.large_sides)
    orientations = CuArray{Int}(cache.mortars.orientations)
    u_upper = CuArray{Float64}(cache.mortars.u_upper)
    u_lower = CuArray{Float64}(cache.mortars.u_lower)
    forward_upper = CuArray{Float64}(dg.mortar.forward_upper)
    forward_lower = CuArray{Float64}(dg.mortar.forward_lower)

    size_arr = CuArray{Float64}(undef, size(u_upper, 2), size(u_upper, 3), size(u_upper, 4))

    prolong_mortars_small2small_kernel =
        @cuda launch = false prolong_mortars_small2small_kernel!(
            u_upper,
            u_lower,
            u,
            neighbor_ids,
            large_sides,
            orientations,
        )
    prolong_mortars_small2small_kernel(
        u_upper,
        u_lower,
        u,
        neighbor_ids,
        large_sides,
        orientations;
        configurator_3d(prolong_mortars_small2small_kernel, size_arr)...,
    )

    prolong_mortars_large2small_kernel =
        @cuda launch = false prolong_mortars_large2small_kernel!(
            u_upper,
            u_lower,
            u,
            forward_upper,
            forward_lower,
            neighbor_ids,
            large_sides,
            orientations,
        )
    prolong_mortars_large2small_kernel(
        u_upper,
        u_lower,
        u,
        forward_upper,
        forward_lower,
        neighbor_ids,
        large_sides,
        orientations;
        configurator_3d(prolong_mortars_large2small_kernel, size_arr)...,
    )

    cache.mortars.u_upper = u_upper # Automatically copy back to CPU
    cache.mortars.u_lower = u_lower # Automatically copy back to CPU

    return nothing
end

# CUDA kernel for calculating surface integrals along axis x, y
function surface_integral_kernel!(du, factor_arr, surface_flux_values)

    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    k = (blockIdx().z - 1) * blockDim().z + threadIdx().z

    if (i <= size(du, 1) && j <= size(du, 2)^2 && k <= size(du, 4))
        j1 = div(j - 1, size(du, 2)) + 1
        j2 = rem(j - 1, size(du, 2)) + 1

        @inbounds begin
            du[i, j1, j2, k] -=
                (
                    surface_flux_values[i, j2, 1, k] * isequal(j1, 1) +
                    surface_flux_values[i, j1, 3, k] * isequal(j2, 1)
                ) * factor_arr[1]
            du[i, j1, j2, k] +=
                (
                    surface_flux_values[i, j2, 2, k] * isequal(j1, size(du, 2)) +
                    surface_flux_values[i, j1, 4, k] * isequal(j2, size(du, 2))
                ) * factor_arr[2]
        end
    end

    return nothing
end

# Launch CUDA kernel to calculate surface integrals
function cuda_surface_integral!(du, mesh::TreeMesh{2}, dg::DGSEM, cache) # surface_integral

    factor_arr = CuArray{Float64}([
        dg.basis.boundary_interpolation[1, 1],
        dg.basis.boundary_interpolation[size(du, 2), 2],
    ])
    surface_flux_values = CuArray{Float64}(cache.elements.surface_flux_values)

    size_arr = CuArray{Float64}(undef, size(du, 1), size(du, 2)^2, size(du, 4))

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

    if (i <= size(du, 1) && j <= size(du, 2)^2 && k <= size(du, 4))
        j1 = div(j - 1, size(du, 2)) + 1
        j2 = rem(j - 1, size(du, 2)) + 1

        @inbounds du[i, j1, j2, k] *= -inverse_jacobian[k]
    end

    return nothing
end

# Launch CUDA kernel to apply Jacobian to reference element
function cuda_jacobian!(du, mesh::TreeMesh{2}, cache)

    inverse_jacobian = CuArray{Float64}(cache.elements.inverse_jacobian)

    size_arr = CuArray{Float64}(undef, size(du, 1), size(du, 2)^2, size(du, 4))

    jacobian_kernel = @cuda launch = false jacobian_kernel!(du, inverse_jacobian)
    jacobian_kernel(du, inverse_jacobian; configurator_3d(jacobian_kernel, size_arr)...)

    return nothing
end

# CUDA kernel for calculating source terms
function source_terms_kernel!(
    du,
    u,
    node_coordinates,
    t,
    equations::AbstractEquations{2},
    source_terms::Function,
)

    j = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    k = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if (j <= size(du, 2)^2 && k <= size(du, 4))
        j1 = div(j - 1, size(du, 2)) + 1
        j2 = rem(j - 1, size(du, 2)) + 1

        u_local = get_nodes_vars(u, equations, j1, j2, k)
        x_local = get_node_coords(node_coordinates, equations, j1, j2, k)

        source_terms_node = source_terms(u_local, x_local, t, equations)

        @inbounds begin
            for ii in axes(du, 1)
                du[ii, j1, j2, k] += source_terms_node[ii]
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
    equations::AbstractEquations{2},
    cache,
)

    return nothing
end

# Launch CUDA kernel to calculate source terms 
function cuda_sources!(du, u, t, source_terms, equations::AbstractEquations{2}, cache)

    node_coordinates = CuArray{Float64}(cache.elements.node_coordinates)

    size_arr = CuArray{Float64}(undef, size(u, 2)^2, size(u, 4))

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

# Pack kernels into `rhs_cpu!()`
#################################################################################
function rhs_cpu!(
    du,
    u,
    t,
    mesh::TreeMesh{2},
    equations,
    initial_condition,
    boundary_conditions,
    source_terms::Source,
    dg::DGSEM,
    cache,
) where {Source}

    reset_du!(du, dg, cache)

    calc_volume_integral!(
        du,
        u,
        mesh,
        have_nonconservative_terms(equations),
        equations,
        dg.volume_integral,
        dg,
        cache,
    )

    prolong2interfaces!(cache, u, mesh, equations, dg.surface_integral, dg)

    calc_interface_flux!(
        cache.elements.surface_flux_values,
        mesh,
        have_nonconservative_terms(equations),
        equations,
        dg.surface_integral,
        dg,
        cache,
    )

    prolong2boundaries!(cache, u, mesh, equations, dg.surface_integral, dg)

    calc_boundary_flux!(
        cache,
        t,
        boundary_conditions,
        mesh,
        equations,
        dg.surface_integral,
        dg,
    )

    prolong2mortars!(cache, u, mesh, equations, dg.mortar, dg.surface_integral, dg)

    calc_mortar_flux!(
        cache.elements.surface_flux_values,
        mesh,
        have_nonconservative_terms(equations),
        equations,
        dg.mortar,
        dg.surface_integral,
        dg,
        cache,
    )

    calc_surface_integral!(du, u, mesh, equations, dg.surface_integral, dg, cache)

    apply_jacobian!(du, mesh, equations, dg, cache)

    calc_sources!(du, u, t, source_terms, equations, dg, cache)

    return nothing
end

function rhs_cpu!(du_ode, u_ode, semi::SemidiscretizationHyperbolic, t)
    @unpack mesh,
    equations,
    initial_condition,
    boundary_conditions,
    source_terms,
    solver,
    cache = semi

    u = wrap_array(u_ode, mesh, equations, solver, cache)
    du = wrap_array(du_ode, mesh, equations, solver, cache)

    rhs_cpu!(
        du,
        u,
        t,
        mesh,
        equations,
        initial_condition,
        boundary_conditions,
        source_terms,
        solver,
        cache,
    )

    return nothing
end

function semidiscretize_cpu(semi::SemidiscretizationHyperbolic, tspan)
    u0_ode = compute_coefficients(first(tspan), semi)

    iip = true
    specialize = SciMLBase.FullSpecialize
    return ODEProblem{iip,specialize}(rhs_cpu!, u0_ode, tspan, semi)
end

# Pack kernels into `rhs_gpu!()`
#################################################################################
function rhs_gpu!(
    du_cpu,
    u_cpu,
    t,
    mesh::TreeMesh{2},
    equations,
    initial_condition,
    boundary_conditions,
    source_terms::Source,
    dg::DGSEM,
    cache,
) where {Source}

    du, u = copy_to_gpu!(du_cpu, u_cpu)

    cuda_volume_integral!(
        du,
        u,
        mesh,
        have_nonconservative_terms(equations),
        equations,
        dg.volume_integral,
        dg,
    )

    cuda_prolong2interfaces!(u, mesh, cache)

    cuda_interface_flux!(mesh, have_nonconservative_terms(equations), equations, dg, cache)

    cuda_prolong2boundaries!(u, mesh, boundary_conditions, cache)

    cuda_boundary_flux!(t, mesh, boundary_conditions, equations, dg, cache)

    cuda_surface_integral!(du, mesh, dg, cache)

    cuda_jacobian!(du, mesh, cache)

    cuda_sources!(du, u, t, source_terms, equations, cache)

    du_computed, _ = copy_to_cpu!(du, u)
    du_cpu .= du_computed

    return nothing
end

function rhs_gpu!(du_ode, u_ode, semi::SemidiscretizationHyperbolic, t)
    @unpack mesh,
    equations,
    initial_condition,
    boundary_conditions,
    source_terms,
    solver,
    cache = semi

    u = wrap_array(u_ode, mesh, equations, solver, cache)
    du = wrap_array(du_ode, mesh, equations, solver, cache)

    rhs_gpu!(
        du,
        u,
        t,
        mesh,
        equations,
        initial_condition,
        boundary_conditions,
        source_terms,
        solver,
        cache,
    )

    return nothing
end

function semidiscretize_gpu(semi::SemidiscretizationHyperbolic, tspan)
    u0_ode = compute_coefficients(first(tspan), semi)

    iip = true
    specialize = SciMLBase.FullSpecialize
    return ODEProblem{iip,specialize}(rhs_gpu!, u0_ode, tspan, semi)
end

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
    equations, solver, cache,)

cuda_prolong2boundaries!(u, mesh,
    boundary_conditions, cache)

cuda_boundary_flux!(t, mesh, boundary_conditions,
    equations, solver, cache)

cuda_prolong2mortars!(u, mesh, solver,
    check_cache_mortars(cache), cache)

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
