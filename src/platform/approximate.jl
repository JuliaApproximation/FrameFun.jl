
function approximate(platform::Platform, param, fun;
            discretizationstyle = DiscretizationStyle(platform),
            solverstyle = SolverStyle(platform),
            verbose = false,
            options...)

    dict = dictionary(platform, param; options...)
    if verbose
        println("Platform approximation:")
        println("- using dictionary $dict")
        println("- using discretization style $discretizationstyle")
        println("- using solver style $solverstyle")
    end
    approximate(discretizationstyle, solverstyle, platform, param, dict, fun; verbose = verbose, options...)
end

approximate(dstyle::DiscretizationStyle, sstyle::SolverStyle, platform::Platform, param, dict, fun; options...) =
    approximate(dstyle, sstyle, dict, fun; options...)

approximate(dict::Dictionary, fun;
            discretizationstyle = default_discretizationstyle(dict),
            solverstyle = default_solverstyle(dict, discretizationstyle),
            options...) =
    approximate(discretizationstyle, solverstyle, dict, fun; options...)


function approximate(dstyle::DiscretizationStyle, sstyle::SolverStyle, dict::Dictionary, fun;
            verbose = false, options...)

    A, B = discretize(dstyle, dict, fun; verbose = verbose, options...)
    C = solve(sstyle, A, B; verbose = verbose, options...)
    DictFun(dict, C)
end

interpolation_grid(dict::Dictionary; options...) = grid(dict)

function discretize(::InterpolationStyle, dict, fun; options...)
    grid = interpolation_grid(dict; options...)
    dicretize_grid(dict, fun, grid)
end

discretize(::GridStyle, dict, fun; grid, options...) = discretize_grid(dict, fun, grid)

function discretize_grid(dict, fun, grid)
    A = evaluation_operator(dict, grid)
    B = sample(grid, fun)
    A, B
end

solve(::InverseStyle, A, B; options...) = inv(A) * B

function directsolver(A; directsolver = :svd, verbose = false, options...)
    if directsolver == :svd
        verbose && println("Using direct solver: SVD")
        SVD_solver(A)
    elseif directsolver == :qr
        verbose && println("Using direct solver: QR")
        QR_solver(A)
    else
        DS = directsolver(A)
        verbose && println("Using direct solver: $DS")
        DS
    end
end

function solve(::DirectStyle, A, B; options...)
    DS = directsolver(A; options...)
    DS * B
end
