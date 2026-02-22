using PlotlyGeometries
using PlotlySupply
using Test

@testset "PlotlyGeometries.jl" begin
    @testset "Base Geometry + Plot Mutation" begin
        c = cuboids([0, 0, 0], [1, 2, 3], "red"; opc=0.5)
        @test c[:type] == "mesh3d"

        gtrans!(c, [1, 2, 3])
        @test isapprox(c.x[1], 0.5; atol=1e-9)
        @test isapprox(c.y[1], 1.0; atol=1e-9)
        @test isapprox(c.z[1], 1.5; atol=1e-9)

        fig = Plot([c], blank_layout())
        add_ref_axes!(fig)
        @test length(fig.data) == 7

        add_arrows!(fig, [0, 0, 0], [1, 0, 0], 1.0, "black")
        add_text!(fig, [0, 0, 0], "origin", "black")
        @test length(fig.data) == 10

        set_view!(fig, 45, 30)
        @test haskey(fig.layout, :scene)
        @test haskey(fig.layout.scene, :camera)
        @test haskey(fig.layout.scene[:camera], :eye)
    end

    @testset "New Primitives" begin
        cyl = cylinders([0, 0, 0], 1.0, 2.0, "z", "royalblue"; tres=24)
        fru = cones([0, 0, 0], 1.0, 0.5, 2.0, "z", "orange"; tres=24)
        con = cones([0, 0, 0], 1.0, 2.0, "z", "orange"; tres=24)
        tor = tori([0, 0, 0], 3.0, 1.0, "z", "teal"; ures=24, vres=12)
        pln = planes([0, 0, 0], [2.0, 3.0], "y", "gray")
        dsk = disks([0, 0, 0], 1.5, "x", "purple"; tres=24)

        for geo in (cyl, fru, con, tor, pln, dsk)
            @test geo[:type] == "mesh3d"
            @test length(geo.i) > 0
            @test length(geo.j) == length(geo.i)
            @test length(geo.k) == length(geo.i)
        end
    end

    @testset "Concave Polygon Triangulation" begin
        pts = [
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
            [1.0, 0.3, 0.0], # concave notch
            [0.0, 1.0, 0.0],
        ]
        poly = polygons(pts, "cyan")
        @test poly[:type] == "mesh3d"
        @test length(poly.i) == length(pts) - 2
    end

    @testset "Group Transforms + Edge Guards" begin
        l1 = lines([0, 0, 0], [1, 0, 0], "red")
        l2 = lines([0, 1, 0], [1, 1, 0], "red")
        geos = [l1, l2]

        gtrans!(geos, [1, 0, 0])
        @test l1.x == [1, 2]
        @test l2.x == [1, 2]

        grot!(geos, 90, [0, 0, 1], [0, 0, 0])
        @test isapprox(l1.x[1], 0.0; atol=1e-9)
        @test isapprox(l1.y[1], 1.0; atol=1e-9)

        @test_throws AssertionError add_arrows!(Plot([cuboids([0, 0, 0], [1, 1, 1])], blank_layout()), [0, 0, 0], [0, 0, 0], 1.0, "black")
        @test_throws AssertionError grot!(cuboids([0, 0, 0], [1, 1, 1]), 30, [0, 0, 0], [0, 0, 0])
        @test_throws AssertionError polygons([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], 3, "blue")
    end

    @testset "STL Import" begin
        # Write a minimal binary STL with 2 triangles (a square)
        stl_path = tempname() * ".stl"
        open(stl_path, "w") do io
            write(io, zeros(UInt8, 80))        # header
            write(io, htol(UInt32(2)))          # 2 triangles
            # triangle 1: (0,0,0) (1,0,0) (1,1,0)
            write(io, htol.(Float32[0, 0, 1]))  # normal
            write(io, htol.(Float32[0, 0, 0]))  # v1
            write(io, htol.(Float32[1, 0, 0]))  # v2
            write(io, htol.(Float32[1, 1, 0]))  # v3
            write(io, UInt16(0))                 # attribute
            # triangle 2: (0,0,0) (1,1,0) (0,1,0)
            write(io, htol.(Float32[0, 0, 1]))
            write(io, htol.(Float32[0, 0, 0]))
            write(io, htol.(Float32[1, 1, 0]))
            write(io, htol.(Float32[0, 1, 0]))
            write(io, UInt16(0))
        end

        m = stlmesh(stl_path, "red")
        @test m[:type] == "mesh3d"
        @test length(m.x) == 6   # 2 triangles * 3 vertices
        @test length(m.i) == 2   # 2 triangles
        @test isapprox(m.x[1], 0.0; atol=1e-6)
        @test isapprox(m.x[2], 1.0; atol=1e-6)
        rm(stl_path)

        # ASCII STL
        stl_ascii_path = tempname() * ".stl"
        open(stl_ascii_path, "w") do io
            println(io, "solid test")
            println(io, "  facet normal 0 0 1")
            println(io, "    outer loop")
            println(io, "      vertex 0.0 0.0 0.0")
            println(io, "      vertex 1.0 0.0 0.0")
            println(io, "      vertex 0.5 1.0 0.0")
            println(io, "    endloop")
            println(io, "  endfacet")
            println(io, "endsolid test")
        end

        ma = stlmesh(stl_ascii_path, "blue")
        @test ma[:type] == "mesh3d"
        @test length(ma.x) == 3
        @test length(ma.i) == 1
        rm(stl_ascii_path)

        @test_throws AssertionError stlmesh("nonexistent.stl")
    end
end
