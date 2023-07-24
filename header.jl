# Import and using packages
using Trixi
using CUDA
using Random, BenchmarkTools, ChangePrecision, OrdinaryDiffEq

using LinearAlgebra: LinearAlgebra, Diagonal, diag, dot, mul!, norm, cross, normalize, I, UniformScaling, det
using Printf: @printf, @sprintf, println
using SparseArrays: AbstractSparseMatrix, AbstractSparseMatrixCSC, sparse, droptol!, rowvals, nzrange, nonzeros, spzeros

using Reexport: @reexport

using MPI: MPI

using SciMLBase: CallbackSet, DiscreteCallback,
    ODEProblem, ODESolution, ODEFunction,
    SplitODEProblem
import SciMLBase: get_du, get_tmp_cache, u_modified!,
    AbstractODEIntegrator, init, step!, check_error,
    get_proposed_dt, set_proposed_dt!,
    terminate!, remake
using CodeTracking: CodeTracking
using ConstructionBase: ConstructionBase
using DiffEqCallbacks: PeriodicCallback, PeriodicCallbackAffect
@reexport using EllipsisNotation
using FillArrays: Ones, Zeros
using ForwardDiff: ForwardDiff
using HDF5: HDF5, h5open, attributes, create_dataset, datatype, dataspace
using IfElse: ifelse
using LinearMaps: LinearMap
using LoopVectorization: LoopVectorization, @turbo, indices
using StaticArrayInterface: static_length
using MuladdMacro: @muladd
using Octavian: Octavian, matmul!
using Polyester: @batch
using OffsetArrays: OffsetArray, OffsetVector
using P4est
using Setfield: @set
using RecipesBase: RecipesBase
using Requires: @require
using Static: Static, One, True, False
@reexport using StaticArrays: SVector
using StaticArrays: StaticArrays, MVector, MArray, SMatrix, @SMatrix
using StrideArrays: PtrArray, StrideArray, StaticInt
@reexport using StructArrays: StructArrays, StructArray
using TimerOutputs: TimerOutputs, @notimeit, TimerOutput, print_timer, reset_timer!
using Triangulate: Triangulate, TriangulateIO, triangulate
export TriangulateIO
using TriplotBase: TriplotBase
using TriplotRecipes: DGTriPseudocolor
@reexport using SimpleUnPack: @unpack
using SimpleUnPack: @pack!

using SummationByPartsOperators: AbstractDerivativeOperator,
    AbstractNonperiodicDerivativeOperator, DerivativeOperator,
    AbstractPeriodicDerivativeOperator, PeriodicDerivativeOperator, grid
import SummationByPartsOperators: integrate, semidiscretize,
    compute_coefficients, compute_coefficients!,
    left_boundary_weight, right_boundary_weight
@reexport using SummationByPartsOperators:
    SummationByPartsOperators, derivative_operator, periodic_derivative_operator,
    upwind_operators

@reexport using StartUpDG: StartUpDG, Polynomial, Gauss, SBP, Line, Tri, Quad, Hex, Tet
using StartUpDG: RefElemData, MeshData, AbstractElemShape

# Include trixi files
include("trixi/src/basic_types.jl")

include("trixi/src/auxiliary/auxiliary.jl")
include("trixi/src/auxiliary/mpi.jl")
include("trixi/src/auxiliary/p4est.jl")
include("trixi/src/equations/equations.jl")
include("trixi/src/meshes/meshes.jl")
include("trixi/src/solvers/solvers.jl")
include("trixi/src/equations/equations_parabolic.jl")
include("trixi/src/semidiscretization/semidiscretization.jl")
include("trixi/src/semidiscretization/semidiscretization_hyperbolic.jl")
include("trixi/src/semidiscretization/semidiscretization_hyperbolic_parabolic.jl")
include("trixi/src/semidiscretization/semidiscretization_euler_acoustics.jl")
include("trixi/src/callbacks_step/callbacks_step.jl")
include("trixi/src/callbacks_stage/callbacks_stage.jl")
include("trixi/src/semidiscretization/semidiscretization_euler_gravity.jl")
include("trixi/src/time_integration/time_integration.jl")

include("trixi/src/auxiliary/special_elixirs.jl")

include("trixi/src/visualization/visualization.jl")