# Remove it after first run to avoid recompilation
#include("header.jl")

# Use the target test header file
include("test/linear_scalar_advection_2d.jl")

# Kernel configurators 
#################################################################################

# CUDA kernel configurator for 1D array computing
function configurator_1d(kernel::CUDA.HostKernel, array::CuArray{Float32,1})
    config = launch_configuration(kernel.fun)

    threads = min(length(array), config.threads)
    blocks = cld(length(array), threads)

    return (threads = threads, blocks = blocks)
end

# CUDA kernel configurator for 2D array computing
function configurator_2d(kernel::CUDA.HostKernel, array::CuArray{Float32,2})
    config = launch_configuration(kernel.fun)

    threads =
        Tuple(fill(Int(floor((min(maximum(size(array)), config.threads))^(1 / 2))), 2))
    blocks = map(cld, size(array), threads)

    return (threads = threads, blocks = blocks)
end

# CUDA kernel configurator for 3D array computing
function configurator_3d(kernel::CUDA.HostKernel, array::CuArray{Float32,3})
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

# CUDA kernel for calculating fluxes along normal direction 1 and 2
function flux_kernel!(
    flux_arr1,
    flux_arr2,
    u,
    equations::AbstractEquations{2},
    flux::Function,
)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    k = (blockIdx().z - 1) * blockDim().z + threadIdx().z

    if (i <= size(u, 1) && j <= size(u, 2)^2 && k <= size(u, 4))
        j1 = div(j - 1, size(u, 2)) + 1
        j2 = rem(j - 1, size(u, 2)) + 1

        u_node = get_nodes_vars(u, equations, j1, j2, k)

        @inbounds begin
            flux_arr1[i, j1, j2, k] = flux(u_node, 1, equations)[i]
            flux_arr2[i, j1, j2, k] = flux(u_node, 2, equations)[i]
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

# Calculate volume integral
function cuda_volume_integral!(
    du,
    u,
    mesh::TreeMesh{2},
    nonconservative_terms,
    equations,
    volume_integral::VolumeIntegralWeakForm,
    dg::DGSEM,
)

    derivative_dhat = CuArray{Float32}(dg.basis.derivative_dhat)
    flux_arr1 = similar(u)
    flux_arr2 = similar(u)
    size_arr = CuArray{Float32}(undef, size(u, 1), size(u, 2)^2, size(u, 4))

    flux_kernel =
        @cuda launch = false flux_kernel!(flux_arr1, flux_arr2, u, equations, flux)
    flux_kernel(
        flux_arr1,
        flux_arr2,
        u,
        equations;
        configurator_3d(flux_kernel, size_arr)...,
    )

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


# Inside `rhs!()` raw implementation
#################################################################################
du, u = copy_to_gpu!(du, u)

cuda_volume_integral!(
    du,
    u,
    mesh,
    have_nonconservative_terms(equations),
    equations,
    solver.volume_integral,
    solver,
)



# For tests
#################################################################################
#= reset_du!(du, solver, cache)

calc_volume_integral!(
    du, u, mesh,
    have_nonconservative_terms(equations), equations,
    solver.volume_integral, solver, cache) =#


#################################################################################

#= nelements(dg::DG, cache) = nelements(cache.elements) =#

#= nelements(elements::ElementContainer2D) = length(elements.cell_ids) =#

#= ntuple(_ -> StaticInt(nnodes(solver)), ndims(mesh))..., nelements(solver, cache) =#

#= unsafe_wrap(Array{eltype(u_ode),ndims(mesh) + 2}, pointer(u_ode),
    (nvariables(equations), ntuple(_ -> nnodes(solver), ndims(mesh))..., nelements(solver, cache))) =#
