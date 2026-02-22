using PlotlySupply
using PlotlyGeometries

# Classic geometries
c1 = cuboids([0, 0, 0], [1, 2, 3], "deepskyblue"; opc=0.28)
e1 = ellipsoids([-6.5, -2.5, 0.0], [2.2, 1.0, 1.4], "purple"; opc=0.22)
s1 = spheres([-2.8, -2.5, 0.0], 1.2, "coral"; opc=0.22)
l1 = lines([-6.5, 2.4, 0.8], [-2.8, -2.5, 0.0], "orange"; style="dash")

# Transform classic geometries
gtrans!(c1, [-6.5, 2.4, 0.8])
grot!(c1, 32, [0, 0, 1])
grot!(e1, [12, 22, 30])

# New primitives
cy1 = cylinders([2.5, 3.0, 0.0], 0.8, 2.5, "z", "steelblue"; opc=0.45, tres=50)
co1 = cones([6.0, 3.0, 0.0], 0.9, 2.0, "z", "tomato"; opc=0.45, tres=50)
fr1 = cones([9.5, 3.0, 0.0], 1.0, 0.45, 2.0, "z", "goldenrod"; opc=0.45, tres=50)
t1 = tori([2.5, -3.0, 0.0], 1.3, 0.45, "z", "teal"; opc=0.35, ures=70, vres=40)
pl1 = planes([6.0, -3.0, 0.0], [2.2, 2.2], "y", "lightgray"; opc=0.40)
dk1 = disks([9.5, -3.0, 0.0], 1.2, "x", "purple"; opc=0.50, tres=60)

# Group transform support on multiple traces
grot!([cy1, co1, fr1, t1, pl1, dk1], [15, 0, 18], [6.0, 0.0, 0.0])

# Create figure with all geometries
fig = plot([c1, e1, s1, l1, cy1, co1, fr1, t1, pl1, dk1], blank_layout())
add_ref_axes!(fig, [0, 0, 0], [2.0, 2.0, 2.0])

# Labels and arrows
add_arrows!(fig, [-6.5, 2.4, 3.2], [0, 0, -1], 0.8, "black")
add_text!(fig, [-6.5, 2.4, 3.8], "cube")

add_arrows!(fig, [-6.0, -0.2, 2.0], [0, -1, -1], 0.8, "black")
add_text!(fig, [-6.0, -0.2, 2.6], "ellipsoid")

add_arrows!(fig, [-2.8, -2.5, 2.0], [0, 0, -1], 0.8, "black")
add_text!(fig, [-2.8, -2.5, 2.6], "sphere")

add_text!(fig, [2.3, 4.8, 0.9], "cylinder", "steelblue")
add_text!(fig, [5.8, 4.8, 0.9], "cone", "tomato")
add_text!(fig, [9.2, 4.8, 0.9], "frustum", "goldenrod")
add_text!(fig, [2.4, -0.8, 0.9], "torus", "teal")
add_text!(fig, [5.9, -0.8, 0.9], "plane", "black")
add_text!(fig, [9.4, -0.8, 0.9], "disk", "purple")

display(fig)
