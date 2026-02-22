module PlotlyGeometries

using PlotlySupply
using Combinatorics
using LinearAlgebra
using BatchAssign

export cuboids, cubes, squares, ellipsoids, spheres, cylinders, cones, tori, planes, disks, lines, polygons, stlmesh,
    grot!, gtrans!, sort_pts, sort_pts!, 
    add_ref_axes!, add_arrows!, add_text!, blank_layout, set_view!

include("api.jl")

end
