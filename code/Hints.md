# Some hints

## Compiling an ITensor System Image

Compile ITensors.jl and make a system image built by the [PackageCompiler.jl](https://github.com/JuliaLang/PackageCompiler.jl) library.

Type the one-line command:

```julia
julia> using ITensors; ITensors.compile()
```

Once ITensors.jl is installed, you can just run this command in an interactive Julia session. It can take several minutes to run, but you only have to run it once for a given version of ITensors.jl. When it is done, it will create a file `sys_itensors.so` in the directory `~/.julia/sysimages/`.

To use the compiled system image together with Julia, run the `julia` command (for interactive mode or scripts) in the following way:

```
$ julia --sysimage ~/.julia/sysimages/sys_itensors.so
```

### Using a Compiled Sysimage in Jupyter or VS Code

To load the ITensor sysimage in VS Code, you can add 
```
"--sysimage ~/.julia/sysimages/sys_itensors.so"
```
as an argument under the `julia.additionalArgs` setting in your Settings.json file.


For more information, see the following [link to the documentation page](https://itensor.github.io/ITensors.jl/dev/getting_started/RunningCodes.html).

## Enabling Debug Checks

ITensor provides some optional checks for common errors, which we call "debug checks".
These can be enabled with the command:
```julia
ITensors.enable_debug_checks()
```
and disabled with the command:
```julia
ITensors.disable_debug_checks()
```

We recommend enabling debug checks when you are developing and testing your code, and then
disabling them when running in production to get the best performance.
