using PlotlySupply
using PlotlyGeometries

# Load tunnel STL
tunnel = stlmesh(joinpath(@__DIR__, "tunnel_curved_arched.stl"), "lightsteelblue"; opc=0.3)

# Floor plane under the tunnel
floor = planes([8, 8, 0], [18, 18], "z", "gainsboro"; opc=0.4)

# -- Items inside the tunnel --

# Barrel near the entrance (y≈2)
barrel = cylinders([15.0, 2.0, 0.4], 0.3, 0.8, "z", "sienna"; opc=0.9, tres=20)

# Crate stack at y≈5
crate1 = cubes([14.1, 5.0, 0.3], 0.6, "goldenrod"; opc=0.9)
crate2 = cubes([14.1, 5.0, 0.9], 0.5, "darkgoldenrod"; opc=0.9)
grot!(crate2, 25, [0, 0, 1])

# Sphere light at y≈8 (hanging from ceiling)
light = spheres([12.6, 8.0, 2.5], 0.15, "gold"; opc=0.95)
wire = lines([12.6, 8.0, 3.0], [12.6, 8.0, 2.65], "dimgray")

# Cone marker at the curve midpoint (y≈11)
cone_marker = cones([10.0, 11.0, 0.0], 0.25, 0.6, "z", "orangered"; opc=0.9, tres=16)

# Torus ring floating near the exit (y≈14)
ring = tori([4.2, 14.5, 1.5], 0.5, 0.12, "x", "mediumseagreen"; opc=0.85, ures=30, vres=16)

# Small sphere at the far end (y≈15)
orb = spheres([5.2, 15.2, 0.4], 0.25, "dodgerblue"; opc=0.85)

# Assemble figure
traces = [tunnel, floor,
          barrel, crate1, crate2,
          light, wire,
          cone_marker, ring, orb]

fig = plot(traces, blank_layout())
add_ref_axes!(fig, [0, 0, 0], [1.5, 1.5, 1.5])

# Labels
add_text!(fig, [15.0, 2.0, 1.5], "barrel", "sienna")
add_text!(fig, [14.1, 5.0, 1.8], "crates", "darkgoldenrod")
add_text!(fig, [12.6, 8.0, 2.0], "light", "goldenrod")
add_text!(fig, [10.0, 11.0, 1.0], "cone", "orangered")
add_text!(fig, [4.2, 14.5, 2.3], "ring", "mediumseagreen")

set_view!(fig, -50, 25)
display(fig)
