# Replicating the `trixi_cuda` Project on Another Machine

Recreating the `trixi_cuda` environment on a different computer is a straightforward procedure based on Julia's `Project.toml` and `Manifest.toml` files. These files capture the precise environment, including all dependencies and their specified versions. Here's a step-by-step guide to do this.

## 1. Transfer Project
Firstly, move your trixi_cuda project to the new machine. This can be accomplished in various ways:
- Using a version control system like Git (by cloning the repository).
- Directly copying the project directory through USB, network transfer, etc.

## 2. Navigate Directory
Once transferred, open a Julia terminal (REPL) on the new machine and navigate to the `trixi_cuda` project's directory,
```Bash
$ cd /path/to/project/trixi_cuda
```

## 3. Activate Project
In the Julia REPL, activate the project using commands like below,
```Julia
julia> using Pkg
julia> Pkg.activate(".")
```

## 4. Instantiate Project
To recreate the environment, use Julia's instantiate function. This function examines the `Project.toml` and `Manifest.toml` files and installs the exact versions of all enumerated packages,
```Julia
julia> Pkg.instantiate()
```
This step may take some time as Julia will fetch and install all required packages.

## 5. Run Project
With the environment now set up, you're ready to run your Julia scripts contained within the `trixi_cuda` project.

By following the above steps, you can guarantee that the Julia environment on the new machine matches the original one. This ensures consistent behavior and performance of the `trixi_cuda` project across various systems.
