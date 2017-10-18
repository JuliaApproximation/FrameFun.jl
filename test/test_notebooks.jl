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
for LINE in eachline(FILE)
    println("Run $(LINE)")
    include(LINE)
    # Following makes things slow but deletes dependencies between notebooks.
    # workspace()
end
close(FILE)
run(`examples/test_notebooks_after.sh`)

println("Create README.md")
run(`jupyter nbconvert --execute --to markdown --output README.md readme.ipynb`)
end
