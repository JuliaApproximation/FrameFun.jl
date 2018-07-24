if VERSION < v"0.7-"
    Pkg.clone("https://github.com/daanhb/BasisFunctions.jl")
    Pkg.checkout("BasisFunctions", "julia-0.7")
    Pkg.build("BasisFunctions")
else
    Pkg.add("https://github.com/daanhb/BasisFunctions.jl#julia-0.7")
end
