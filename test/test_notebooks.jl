module test_nodebooks

using Conda
using DomainSets
using BasisFunctions
using FrameFun
using Plots

function delimit(s::AbstractString)
    println()
    println("############")
    println("# ",s)
    println("############")
end





jupyter = Sys.which("jupyter")
jupyter = nothing
if jupyter == nothing
    jupyter = Sys.which(Conda.PYTHONDIR * "/jupyter")
    if jupyter == nothing
        try
            @info "Trying to install jupyter using the Conda package"
            Conda.add("jupyter")
            @show Conda.PYTHONDIR
            @show Base.Filesystem.readdir(Conda.PYTHONDIR)
            jupyter = Sys.which(Conda.PYTHONDIR * "/jupyter")
        catch error
            @warn Conda could not add jupyter
            bt = catch_backtrace()
            msg = sprint(showerror, error, bt)
            @info msg
        end
    end
end
if jupyter == nothing
    @warn "Jupyter not installed and no internet connection to install, abandoning notebook tests."
else
    FRAMEFUNSRC = pathof(FrameFun)
    FRAMEFUNPATH = splitdir(splitdir(FRAMEFUNSRC)[1])[1]
    delimit("Notebooks")

    run(`$jupyter nbconvert --to script $FRAMEFUNPATH'/examples/*.ipynb'`)
    a = readdir("$FRAMEFUNPATH/examples/")
    scripts = a[occursin.(".jl",a)]
    try
        for LINE in scripts
            println("Run $(LINE)")
            include("$FRAMEFUNPATH/examples/"*LINE)
        end
    finally
        for LINE in scripts
            rm("$FRAMEFUNPATH/examples/"*LINE)
        end
    end
end
end
