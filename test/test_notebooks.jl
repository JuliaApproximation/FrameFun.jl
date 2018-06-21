module test_nodebooks
using Domains
using BasisFunctions
using FrameFun
using Plots

function delimit(s::AbstractString)
    println()
    println("############")
    println("# ",s)
    println("############")
end

delimit("Notebooks")
run(`examples/test_notebooks.sh`)

FILE = open("notebookscripts")
try
    for LINE in eachline(FILE)
        println("Run $(LINE)")
        include(LINE)
        # Following makes things slow but deletes dependencies between notebooks.
        # workspace()
    end
catch y
    if isa(y, OutOfMemoryError)
        travis==false
        try
            travis = (ENV["TRAVIS"]=="TRUE")
        end
        if travis
            warn("Our of memory")
        else rethrow(y)
        end
    else
        rethrow(y)
    end
finally
    close(FILE)
    run(`examples/test_notebooks_after.sh`)
end

end
