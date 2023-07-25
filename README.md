# linear_advection_cuda

- The folder `trixi` stores folders from `Trixi.jl`, specifically `trixi/src`, `trixi/examples`, and `trixi/test`.
- The file `header.jl` can be used as test environment initializer for running tests.
- The folder `test` contains the header part of tests for different equations.
- The file `cuda_dg_1d.jl`, `cuda_dg_2d.jl`, and `cuda_dg_3d.jl` run tests for prototyping 1D, 2D, and 3D GPU code for `rhs!()` functions.
- The file `simple_kernels.jl` creates sample kernels for running 1D, 2D, and 3D GPU code.
- The folder `docs` contains useful resources for this project.
- The folder `profile` contains contents about how to profile GPU kernels.

TODO List:
- 1D: 1) `nonconservative_terms::True`, 2) `calc_volume_integral!()` - `volume_integral::VolumeIntegralShockCapturingHG`
- 2D: 2) `nonconservative_terms::True`, 2) `calc_volume_integral!()` - `volume_integral::VolumeIntegralShockCapturingHG`, 3) `prolong2mortars!()` and `calc_mortar_flux!()`
- 3D: 1) `nonconservative_terms::True`, 2) `calc_volume_integral!()` - `volume_integral::VolumeIntegralShockCapturingHG`, 3) `prolong2mortars!()` and `calc_mortar_flux!()` and `calc_boundary_flux!()` (issue #6)
- Run something people coming from numerical mathematics care about
- Setup some regression tests comparing your results to CPU results obtained with the current version of Trixi.jl
- Create meaningful benchmarks that can be used to tune the implementations later