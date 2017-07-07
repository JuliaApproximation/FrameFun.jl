# A test suite for support functions
module test_suite_support

using Domains
using BasisFunctions
using FrameFun
using Base.Test

## Settings

# Test fourier extensions for all parameters
Extensive = false
########
# Auxiliary functions
########



function delimit(s::AbstractString)
    println("############")
    println("# ",s)
    println("############")
end

#######
# Testing
#######


function test_subgrids()
    delimit("Grid functionality")
    n = 20
    grid1 = EquispacedGrid(n, -1.0, 1.0)
    subgrid1 = MaskedGrid(grid1, interval(-0.5, 0.7))
    subgrid2 = IndexSubGrid(grid1, 4:12)

    G1 = EquispacedGrid(n, -1.0, 1.0)
    G2 = EquispacedGrid(n, -1.0, 1.0)
    TensorG = G1 ⊗ G2

    C = disk(1.0)
    circle_grid = MaskedGrid(TensorG, C)
    @testset begin

        @test (length(circle_grid)/length(TensorG)-pi*0.25) < 0.01

        G1s = IndexSubGrid(G1,2:4)
        G2s = IndexSubGrid(G2,3:5)
        TensorGs = G1s ⊗ G2s
        @test G1s[1] == G1[2]
        @test G2s[1] == G2[3]
        @test TensorGs[1,1] == [G1[2],G2[3]]
    end

    # Generic tests for the subgrids
    @testset "result" for (grid,subgrid) in ( (grid1,subgrid1), (grid1,subgrid2), (TensorG, circle_grid))
        # print("Subgrid is ")
        # println(typeof(subgrid))
        # Count the number of elements in the subgrid
        cnt = 0
        for i in 1:length(grid)
            if is_subindex(i, subgrid)
                cnt += 1
            end
        end
        @test cnt == length(subgrid)

        space = DiscreteGridSpace(grid)
        subspace = DiscreteGridSpace(subgrid)
        R = restriction_operator(space, subspace)
        E = extension_operator(subspace, space)

        e = random_expansion(subspace)
        e_ext = E * e
        # Are the elements in the right place?
        cnt = 0
        diff = 0.0
        for i in 1:length(grid)
            if is_subindex(i, subgrid)
                cnt += 1
                diff += abs(e[cnt] - e_ext[i])
            else
                # all other entries should be zero
                diff += abs(e_ext[i])
            end
        end
        @test diff < 1e-6

        e_rest = R * e_ext
        @test sum([abs(e[i]-e_rest[i]) for i in 1:length(e)]) < 1e-6
    end
end

function test_randomgrids()
    println("Random grids:")
    @testset begin
        d = disk()
        g = randomgrid(d, 10)
        @test length(g) == 10
        @test length(eltype(g)) == ndims(d)
        @test reduce(&, [x ∈ d for x in g])

        g2 = randomgrid(disk(BigFloat), 10)
        @test eltype(g2[1]) == BigFloat

        g3 = randomgrid(interval(0.0, 1.0), 10)
        @test length(g3) == 10
        # 1D is a special case where we don't use vectors of length 1
        @test eltype(g3) == Float64

        box1 = rectangle(0.0, 1.0, 0.0, 1.0)
        p1 = randompoint(box1)
        @test length(p1) == ndims(box1)
        @test p1 ∈ box1
        box2 = interval(0.0, 1.0)
        p2 = randompoint(box2)
        @test typeof(p2) == Float64
        @test p2 ∈ box2
    end
end



test_subgrids()
test_randomgrids()


end
