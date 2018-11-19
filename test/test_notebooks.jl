module test_nodebooks

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


jupyter = homedir()*"/.julia/packages/Conda/S0nV/deps/usr/bin/jupyter"



delimit("Notebooks")
try
    cd("test/")
catch y
    nothing
finally
    run(`$jupyter nbconvert --to script '../examples/*.ipynb'`)
    run(`mkdir -p scripts/`)
    run(`./findscripts.sh`)
    cd("..")
end

FILE = open("notebookscripts")
try
    for LINE in eachline(FILE)
        println("Run $(LINE)")
        include(LINE)
        # Following makes things slow but deletes dependencies between notebooks.
        # workspace()
    end
catch y
    rethrow(y)
finally
    close(FILE)
    run(`examples/test_notebooks_after.sh`)
end

end
