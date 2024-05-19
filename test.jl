using CUDA, Test, BenchmarkTools, StaticArrays


a = CUDA.fill(SVector(1.0, 2.0, 3.0), (2, 2))
b = CuArray(fill(SVector(1.0, 2.0, 3.0), (2, 2)))
c = CUDA.fill(zeros(3), (2, 2))
d = CuArray(fill(zeros(3), (2, 2)))
