# A test suite for support functions
module test_suite_support


using BasisFunctions
using FrameFun
using Base.Test
FE = FrameFun
BA = BasisFunctions

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
    subgrid1 = FE.MaskedGrid(grid1, Interval(-0.5, 0.7))
    subgrid2 = FE.IndexSubGrid(grid1, 4, 12)

    G1 = EquispacedGrid(n, -1.0, 1.0)
    G2 = EquispacedGrid(n, -1.0, 1.0)
    TensorG = G1 ⊗ G2

    C = Disk(1.0)
    circle_grid = FE.MaskedGrid(TensorG, C)
    @testset begin

        @test (length(circle_grid)/length(TensorG)-pi*0.25) < 0.01

        G1s = FE.IndexSubGrid(G1,2,4)
        G2s = FE.IndexSubGrid(G2,3,5)
        TensorGs = G1s ⊗ G2s
        @test G1s[1] == G1[2]
        @test G2s[1] == G2[3]
        @test TensorGs[1,1] == [G1[2],G2[3]]
    end

    # Generic tests for the subgrids
    @testset "result" for (grid,subgrid) in ( (grid1,subgrid1), (grid1,subgrid2), (TensorG, circle_grid))
        print("Subgrid is ")
        println(typeof(subgrid))
        # Count the number of elements in the subgrid
        cnt = 0
        for i in 1:length(grid)
            if i ∈ subgrid
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
            if i ∈ subgrid
                cnt += 1
                diff += abs(e[cnt] - e_ext[i])
            end
        end
        @test diff < 1e-6

        e_rest = R * e_ext
        @test sum([abs(e[i]-e_rest[i]) for i in 1:length(e)]) < 1e-6
    end
end





test_subgrids()


delimit("Domains")

@testset begin
    # Interval
    Intervala=Interval(-1.0,1.0)
    Intervala=Intervala+1
    Intervala=1+Intervala
    @test FE.left(Intervala)==1
    @test FE.left(2*Intervala)==2
    @test FE.right(Intervala/4)==0.75
    # Disk
    C = Disk(2.0)
    @test FE.in([1.4, 1.4], C)
    @test !FE.in([1.5, 1.5], C)
    @test FE.boundingbox(C) == FE.BBox((-2.0,-2.0),(2.0,2.0))
    # This is certainly unwanted behavior! Due to method inheritance
    @test typeof(1.2*C)==typeof(C*1.2)
    # This is due to a wrong implementation in Scaled Domain
    @test FE.in([1.5,1.5],1.2*C)
    @test FE.in([1.5,1.5],C*1.2)

    #Square
    D = Cube(2)
    @test FE.in([0.9, 0.9],D)
    @test !FE.in([1.1, 1.1],D)
    @test FE.boundingbox(D)==FE.BBox((-1,-1),(1,1))
    DS=FE.union(D,C)
    #Cube
    D=Cube((-1.5,0.5,-3.0),(2.2,0.7,-1.0))
    @test FE.in([0.9, 0.6, -2.5],D)
    @test !FE.in([0.0, 0.6, 0.0],D)
    @test FE.boundingbox(D)==FE.BBox((-1.5,0.5,-3.0),(2.2,0.7,-1.0))

    #Ball
    S=Ball(2.0)
    @test FE.in([1.9,0.0,0.0],S)
    @test FE.in([0,-1.9,0.0],S)
    @test FE.in([0.0,0.0,-1.9],S)
    @test !FE.in([1.9,1.9,0.0],S)
    @test FE.boundingbox(S)==FE.BBox((-2.0,-2.0,-2.0),(2.0,2.0,2.0))
    # joint domain
    DS=FE.union(D,S)
    @test FE.in([0.0,0.6,0.0],DS)
    @test FE.in([0.9, 0.6, -2.5],DS)
    @test FE.boundingbox(DS)==FE.BBox((-2.0,-2.0,-3.0),(2.2,2.0,2.0))
    # domain intersection
    DS=FE.intersect(D,S)
    @test !FE.in([0.0,0.6,0.0],DS)
    @test FE.in([0.2, 0.6, -1.1],DS)
    @test FE.boundingbox(DS)==FE.BBox((-1.5,0.5,-2.0),(2.0,0.7,-1.0))
    # domain difference
    DS=D-S
    @test FE.boundingbox(DS)==FE.boundingbox(D)
    # TensorProductDomain 1
    T=FE.tensorproduct(Interval(-1.0,1.0),2)
    FE.in([0.5,0.5],T)
    @test FE.in([0.5,0.5],T)
    @test !FE.in([-1.1,0.3],T)

    # TensorProductDomain 2
    T=FE.TensorProductDomain(Disk(1.05),Interval(-1.0,1.0))
    @test FE.in([0.5,0.5,0.8],T)
    @test !FE.in([-1.1,0.3,0.1],T)

end

end
