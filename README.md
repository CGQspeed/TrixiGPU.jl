# Provide GPU Support for Trixi.jl via CUDA
This project is undertaken as part of the [Google Summer of Code 2023](https://summerofcode.withgoogle.com/) program and is in the developing and testing phase. Please check the [Project Summary](https://gist.github.com/huiyuxie/44b561f9f854aada98fdb37036081454) for future steps. 

## GPU Strategy Overview
At present, two primary strategies are being explored for leveraging GPU capabilities:

**Utilizing CUDA.jl (Julia):** This approach employs CUDA.jl, which is a Julia package providing an interface to NVIDIA's CUDA APIs. CUDA.jl is effectively a high-level abstraction over the CUDA toolkit, allowing developers to write GPU code using Julia's intuitive syntax. It abstracts away many of the complexities of direct CUDA programming and makes it easier to integrate GPU acceleration into Julia programs. By using CUDA.jl, we're able to directly tap into the power of the GPU without leaving the Julia environment.

**Direct CUDA through C++:** This strategy involves writing GPU code using the native CUDA toolkit in C++. The advantage of this method is the potential for finer control over GPU operations, possibly leading to optimizations that might be more challenging to achieve with a higher-level API. To bridge the gap between Julia and C++, Julia's C and C++ FFI (Foreign Function Interface) is employed. This interface allows data and structures to be transferred seamlessly between Julia and C++. Hence, with this approach, GPU operations are coded and executed in C++, with Julia serving as the orchestrating layer, managing data exchange and other high-level tasks.  

## Project Directory Structure
- The folder `trixi` stores folders from `Trixi.jl`, specifically `trixi/src`, `trixi/examples`, and `trixi/test`.
- The folder `docs` contains useful resources and instructions for this project.
- The folder `profile` contains contents about how to profile GPU kernels.
- The folder `cuda_julia` contains primary files (`cuda_dg_1d.jl`, `cuda_dg_2d.jl`, and `cuda_dg_3d.jl`) that are used to test the performance of `rhs!()` functions implemented in CUDA.jl.
- The folder `cuda_cpp` contains primary files (`cuda_dg_1d.cu`, `cuda_dg_2d.cu`, and `cuda_dg_3d.cu`) that are used to test the performance of `rhs!()` functions implemented in CUDA/C++.

## Kernels to be Implemented
- 1D Kernels: 1) `calc_volume_integral!()` - `volume_integral::VolumeIntegralShockCapturingHG`
- 2D Kernels: 1) `calc_volume_integral!()` - `volume_integral::VolumeIntegralShockCapturingHG`, 2) `calc_mortar_flux!()`
- 3D Kernels: 1) `calc_volume_integral!()` - `volume_integral::VolumeIntegralShockCapturingHG`, 2) `prolong2mortars!()`, 3) `calc_mortar_flux!()` 

## How to Show Your Support
If you found this project interesting and inspiring, kindly give it a star. Your support means a lot to us! ❤️
