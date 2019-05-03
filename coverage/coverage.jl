using Pkg
Pkg.instantiate()

# Only run coverage from linux nightly build on travis.
get(ENV, "TRAVIS_OS_NAME", "")       == "linux"   || exit()
get(ENV, "TRAVIS_JULIA_VERSION", "") == "nightly" || exit()

using Coverage

cd(joinpath(dirname(@__FILE__), "..")) do
    Coveralls.submit(process_folder())
    Codecov.submit(process_folder())
end
