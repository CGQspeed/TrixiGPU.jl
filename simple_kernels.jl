### This is a test file for running CUDA simple kernels, which use different combinations 
### of 1, 2, and 3 dimensions for both blocks and grids to conduct computing on 1D, 2D, and 
### 3D arrays.
### Note that not all the combinations are implemented here due to the similarity.

using CUDA, Test

#--------------------------------------------------------------------- Kernels for 1D array
# In total there are 3 choices for 1D array computing: 
# 1) 1D block; 2) 1D grid; 2) 1D block, 1D grid

# 1D block 
function gpu_add1!(y, x)
    index = threadIdx().x
    stride = blockDim().x

    for i ∈ index:stride:length(y)
        @inbounds y[i] += x[i]
    end

    return nothing
end

N = 2^20
x = CUDA.fill(1.0f0, N)
y = CUDA.fill(2.0f0, N)
@cuda threads = 256 gpu_add1!(y, x)
@test all(Array(y) .== 3.0f0)

# 1D grid
function gpu_add2!(y, x)
    index = blockIdx().x
    stride = gridDim().x

    for i ∈ index:stride:length(y)
        @inbounds y[i] += x[i]
    end

    return nothing
end

N = 2^20
x = CUDA.fill(1.0f0, N)
y = CUDA.fill(2.0f0, N)
@cuda blocks = 256 gpu_add2!(y, x)
@test all(Array(y) .== 3.0f0)

# 1D block, 1D grid
function gpu_add3!(y, x)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = gridDim().x * blockDim().x

    for i ∈ index:stride:length(y)
        @inbounds y[i] += x[i]
    end

    return nothing
end

N = 2^20
x = CUDA.fill(1.0f0, N)
y = CUDA.fill(2.0f0, N)
numblocks = ceil(Int, N / 256)
@cuda threads = 256 blocks = numblocks gpu_add3!(y, x)
@test all(Array(y) .== 3.0f0)

#--------------------------------------------------------------------- Kernels for 2D array
# In total there are 5 choices for 2D array computing:
# 1) 2D block; 2) 2D grid; 3) 2D block, 1D grid; 4) 1D block, 2D grid; 5) 2D block, 2D grid

# 2D block
function gpu_add4!(y, x)
    index_i = threadIdx().x
    index_j = threadIdx().y

    stride_i = blockDim().x
    stride_j = blockDim().y

    for i ∈ index_i:stride_i:size(y, 1)
        for j ∈ index_j:stride_j:size(y, 2)
            @inbounds y[i, j] += x[i, j]
        end
    end

    return nothing
end

N = 2^10
x = CUDA.ones(N, N)
y = CUDA.zeros(N, N)
@cuda threads = (16, 16) gpu_add4!(y, x)
@test all(Array(y) .== 1.0f0)

# 2D block, 1D grid
function gpu_add5!(y, x)
    index_i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    index_j = threadIdx().y

    stride_i = gridDim().x * blockDim().x
    stride_j = blockDim().y

    for i ∈ index_i:stride_i:size(y, 1)
        for j ∈ index_j:stride_j:size(y, 2)
            @inbounds y[i, j] += x[i, j]
        end
    end

    return nothing
end

N = 2^10
x = CUDA.ones(N, N)
y = CUDA.zeros(N, N)
numblocks = ceil(Int, N / 256)
@cuda threads = (16, 16) blocks = numblocks gpu_add5!(y, x)
@test all(Array(y) .== 1.0f0)

# 2D block, 2D grid
function gpu_add6!(y, x)
    index_i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    index_j = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    stride_i = gridDim().x * blockDim().x
    stride_j = gridDim().y * blockDim().y

    for i ∈ index_i:stride_i:size(y, 1)
        for j ∈ index_j:stride_j:size(y, 2)
            @inbounds y[i, j] += x[i, j]
        end
    end

    return nothing
end

N = 2^10
x = CUDA.ones(N, N)
y = CUDA.zeros(N, N)
numblocks = ceil(Int, N / 256)
@cuda threads = (16, 16) blocks = (numblocks, numblocks) gpu_add6!(y, x)
@test all(Array(y) .== 1.0f0)

#--------------------------------------------------------------------- Kernels for 3D array
# In total there are 7 choices for 3D array computing:
# 1) 3D block; 2) 3D grid; 3) 3D block, 1D grid; 4) 1D block, 3D grid; 5) 3D block, 2D grid;
# 6) 2D block, 3D grid; 7) 3D block, 3D grid

# 3D block
function gpu_add7!(y, x)
    index_i = threadIdx().x
    index_j = threadIdx().y
    index_k = threadIdx().z

    stride_i = blockDim().x
    stride_j = blockDim().y
    stride_k = blockDim().z

    for i ∈ index_i:stride_i:size(y, 1)
        for j ∈ index_j:stride_j:size(y, 2)
            for k ∈ index_k:stride_k:size(y, 3)
                @inbounds y[i, j, k] += x[i, j, k]
            end
        end
    end

    return nothing
end

N = 2^8
x = CUDA.ones(N, N, N)
y = CUDA.zeros(N, N, N)
@cuda threads = (8, 8, 8) gpu_add7!(y, x)
@test all(Array(y) .== 1.0f0)

# 3D block, 1D grid
function gpu_add8!(y, x)
    index_i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    index_j = threadIdx().y
    index_k = threadIdx().z

    stride_i = gridDim().x * blockDim().x
    stride_j = blockDim().y
    stride_k = blockDim().z

    for i ∈ index_i:stride_i:size(y, 1)
        for j ∈ index_j:stride_j:size(y, 2)
            for k ∈ index_k:stride_k:size(y, 3)
                @inbounds y[i, j, k] += x[i, j, k]
            end
        end
    end

    return nothing
end

N = 2^8
x = CUDA.ones(N, N, N)
y = CUDA.zeros(N, N, N)
@cuda threads = (8, 8, 8) blocks = 4 gpu_add8!(y, x)
@test all(Array(y) .== 1.0f0)

# 3D block, 2D grid
function gpu_add9!(y, x)
    index_i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    index_j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    index_k = threadIdx().z

    stride_i = gridDim().x * blockDim().x
    stride_j = gridDim().y * blockDim().y
    stride_k = blockDim().z

    for i ∈ index_i:stride_i:size(y, 1)
        for j ∈ index_j:stride_j:size(y, 2)
            for k ∈ index_k:stride_k:size(y, 3)
                @inbounds y[i, j, k] += x[i, j, k]
            end
        end
    end

    return nothing
end

N = 2^8
x = CUDA.ones(N, N, N)
y = CUDA.zeros(N, N, N)
@cuda threads = (8, 8, 8) blocks = (4, 4) gpu_add9!(y, x)
@test all(Array(y) .== 1.0f0)

# 3D block, 3D grid
function gpu_add10!(y, x)
    index_i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    index_j = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    index_k = (blockIdx().z - 1) * blockDim().z + threadIdx().z

    stride_i = gridDim().x * blockDim().x
    stride_j = gridDim().y * blockDim().y
    stride_k = gridDim().z * blockDim().z

    for i ∈ index_i:stride_i:size(y, 1)
        for j ∈ index_j:stride_j:size(y, 2)
            for k ∈ index_k:stride_k:size(y, 3)
                @inbounds y[i, j, k] += x[i, j, k]
            end
        end
    end

    return nothing
end

N = 2^8
x = CUDA.ones(N, N, N)
y = CUDA.zeros(N, N, N)
@cuda threads = (8, 8, 8) blocks = (4, 4, 4) gpu_add10!(y, x)
@test all(Array(y) .== 1.0f0)
