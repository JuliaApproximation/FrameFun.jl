module test_nodebooks

using Conda: PYTHONDIR
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

jupyter = Conda.PYTHONDIR * "jupyter"

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
