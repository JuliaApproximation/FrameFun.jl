
ENV["JULIA_DEBUG"] = "all"

include("test_suite.jl")
include("test_scenariolist.jl")
include("test_suite_applications.jl")
include("test_suite_domains.jl")
include("test_platforms.jl")
include("test_plots.jl")

include("test_notebooks.jl")
# include("create_readme.jl")


println("All tests succeeded!!!")
