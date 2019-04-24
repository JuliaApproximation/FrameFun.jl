
run(`rm -f '~/.jupyter/jupyter_nbconvert_config.json'`)
using Pkg
Pkg.activate(splitdir(@__DIR__)[1])
pkg"add https://github.com/daanhb/BasisFunctions.jl"
