const PlotContainer = Union{Plot, SyncPlot}
const _GEOM_EPS = 1e-10

_cross2d(ax::Real, ay::Real, bx::Real, by::Real) = ax * by - ay * bx

function _random_rgb()
    @all r g b = round(Int, rand() * 255)
    return "rgb($r, $g, $b)"
end

function _validate_axis(mode::String)
    @assert mode in ("x", "y", "z")
    return mode
end

function _orient_point(x::Real, y::Real, z::Real, mode::String)
    if mode == "z"
        return x, y, z
    elseif mode == "x"
        # Rotate around y by +90 deg so z-axis geometry points along x.
        return z, y, -x
    else
        # Rotate around x by -90 deg so z-axis geometry points along y.
        return x, z, -y
    end
end

function _orient_and_shift!(x::Vector{<:Real}, y::Vector{<:Real}, z::Vector{<:Real}, origin::Vector{<:Real}, mode::String)
    @assert length(origin) == 3
    _validate_axis(mode)

    @inbounds for n in eachindex(x)
        xo, yo, zo = _orient_point(x[n], y[n], z[n], mode)
        x[n] = xo + origin[1]
        y[n] = yo + origin[2]
        z[n] = zo + origin[3]
    end
    return nothing
end

function _validate_points3d(pts::Vector)
    @assert all(length.(pts) .== 3)
    for vec in pts
        for num in vec
            @assert isreal(num)
        end
    end
    return nothing
end

function _project_polygon_to_2d(pts::Vector)
    N = length(pts)
    @assert N >= 3

    normal = zeros(Float64, 3)
    found_normal = false
    for a = 1:N-2
        for b = a+1:N-1
            for c = b+1:N
                v1 = pts[b] .- pts[a]
                v2 = pts[c] .- pts[a]
                cand = cross(v1, v2)
                if norm(cand) > _GEOM_EPS
                    normal .= cand ./ norm(cand)
                    found_normal = true
                    break
                end
            end
            found_normal && break
        end
        found_normal && break
    end
    @assert found_normal "Polygon points are collinear; triangulation is undefined."

    origin = pts[1]
    e1 = nothing
    for n = 2:N
        vec = pts[n] .- origin
        vec_plane = vec .- dot(vec, normal) * normal
        if norm(vec_plane) > _GEOM_EPS
            e1 = vec_plane ./ norm(vec_plane)
            break
        end
    end
    @assert !isnothing(e1) "Failed to construct a polygon plane basis."
    e2 = cross(normal, e1)

    uv = Vector{Vector{Float64}}(undef, N)
    for n in eachindex(pts)
        vec = pts[n] .- origin
        uv[n] = [dot(vec, e1), dot(vec, e2)]
    end
    return uv
end

function _polygon_signed_area2d(uv::Vector{<:Vector{<:Real}})
    area = 0.0
    N = length(uv)
    for n in 1:N
        m = n == N ? 1 : n + 1
        area += uv[n][1] * uv[m][2] - uv[m][1] * uv[n][2]
    end
    return area / 2
end

function _point_in_triangle2d(p::Vector{<:Real}, a::Vector{<:Real}, b::Vector{<:Real}, c::Vector{<:Real})
    c1 = _cross2d(b[1] - a[1], b[2] - a[2], p[1] - a[1], p[2] - a[2])
    c2 = _cross2d(c[1] - b[1], c[2] - b[2], p[1] - b[1], p[2] - b[2])
    c3 = _cross2d(a[1] - c[1], a[2] - c[2], p[1] - c[1], p[2] - c[2])
    has_neg = c1 < -_GEOM_EPS || c2 < -_GEOM_EPS || c3 < -_GEOM_EPS
    has_pos = c1 > _GEOM_EPS || c2 > _GEOM_EPS || c3 > _GEOM_EPS
    return !(has_neg && has_pos)
end

function _triangulate_polygon_indices(pts::Vector)
    @assert length(pts) >= 3
    uv = _project_polygon_to_2d(pts)
    idx = collect(1:length(pts))
    if _polygon_signed_area2d(uv) < 0
        reverse!(idx)
    end

    tris = Vector{NTuple{3, Int}}()
    iter_guard = 0
    while length(idx) > 3
        iter_guard += 1
        @assert iter_guard <= 10 * length(pts) "Failed to triangulate polygon. Ensure points form a simple coplanar boundary."

        ear_found = false
        M = length(idx)
        for m in 1:M
            ia = idx[m == 1 ? M : m - 1]
            ib = idx[m]
            ic = idx[m == M ? 1 : m + 1]

            a = uv[ia]
            b = uv[ib]
            c = uv[ic]

            convex = _cross2d(b[1] - a[1], b[2] - a[2], c[1] - a[1], c[2] - a[2]) > _GEOM_EPS
            convex || continue

            any_inside = false
            for id in idx
                if id == ia || id == ib || id == ic
                    continue
                end
                if _point_in_triangle2d(uv[id], a, b, c)
                    any_inside = true
                    break
                end
            end
            any_inside && continue

            push!(tris, (ia, ib, ic))
            deleteat!(idx, m)
            ear_found = true
            break
        end

        @assert ear_found "Failed to triangulate polygon. Ensure points are ordered around a non-self-intersecting boundary."
    end

    push!(tris, (idx[1], idx[2], idx[3]))
    return tris
end

function _append_polygon_mesh!(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64},
        i::Vector{Int}, j::Vector{Int}, k::Vector{Int}, pts::Vector, offset::Int)
    @assert length(pts) >= 3
    pts_ordered = sort_pts(pts)
    tris = _triangulate_polygon_indices(pts_ordered)

    for pt in pts_ordered
        push!(x, pt[1])
        push!(y, pt[2])
        push!(z, pt[3])
    end

    for tri in tris
        push!(i, offset + tri[1] - 1)
        push!(j, offset + tri[2] - 1)
        push!(k, offset + tri[3] - 1)
    end
    return offset + length(pts_ordered)
end

function _is_stl_binary(filepath::String)
    filesize(filepath) < 84 && return false
    open(filepath) do io
        read(io, 80)
        ntri = ltoh(read(io, UInt32))
        return filesize(filepath) == 84 + ntri * 50
    end
end

function _parse_stl_binary(io::IO)
    read(io, 80)
    ntri = ltoh(read(io, UInt32))
    x = Float64[]; y = Float64[]; z = Float64[]
    i = Int[]; j = Int[]; k = Int[]
    for t in 1:ntri
        read(io, 12)  # skip normal
        for v in 1:3
            push!(x, Float64(ltoh(read(io, Float32))))
            push!(y, Float64(ltoh(read(io, Float32))))
            push!(z, Float64(ltoh(read(io, Float32))))
        end
        idx = (t - 1) * 3
        push!(i, idx); push!(j, idx + 1); push!(k, idx + 2)
        read(io, 2)  # skip attribute
    end
    return x, y, z, i, j, k
end

function _parse_stl_ascii(io::IO)
    x = Float64[]; y = Float64[]; z = Float64[]
    i = Int[]; j = Int[]; k = Int[]
    tri_count = 0
    for line in eachline(io)
        tokens = split(strip(line))
        if length(tokens) == 4 && tokens[1] == "vertex"
            push!(x, parse(Float64, tokens[2]))
            push!(y, parse(Float64, tokens[3]))
            push!(z, parse(Float64, tokens[4]))
            if mod(length(x), 3) == 0
                idx = tri_count * 3
                push!(i, idx); push!(j, idx + 1); push!(k, idx + 2)
                tri_count += 1
            end
        end
    end
    return x, y, z, i, j, k
end

"""
    cuboids(origin::Vector{<:Real}, dimension::Vector{<:Real}, color::String=""; opc::Real=1)

    Creates a 3D box mesh centered at the given origin with specified dimensions and color.

    # Arguments
    - `origin::Vector{<:Real}`: A vector of three Reals specifying the center of the box.
    - `dimension::Vector{<:Real}`: A vector of three Reals specifying the dimensions (width, height, depth) of the box.
    - `color::String`: A string specifying the color of the box.
    
    # Keywords
    - `opc`: (optional) A Real specifying the opacity of the box. Default is 1.                                   
"""
function cuboids(origin::Vector{<:Real}, dimension::Vector{<:Real}, color::String=""; opc::Real=1)
    @assert length(origin) == 3
    @assert length(dimension) == 3

    if color == ""
        @all r g b = round(Int, rand() * 255)
        color = "rgb($r, $g, $b)"
    end

    x1 = origin[1] - dimension[1] / 2
    x2 = origin[1] + dimension[1] / 2
    y1 = origin[2] - dimension[2] / 2
    y2 = origin[2] + dimension[2] / 2
    z1 = origin[3] - dimension[3] / 2
    z2 = origin[3] + dimension[3] / 2

    x = [x1, x1, x2, x2, x1, x1, x2, x2]
    y = [y1, y2, y2, y1, y1, y2, y2, y1]
    z = [z1, z1, z1, z1, z2, z2, z2, z2]
    i = [7, 0, 0, 0, 4, 4, 6, 6, 4, 0, 3, 2]
    j = [3, 4, 1, 2, 5, 6, 5, 2, 0, 1, 6, 3]
    k = [0, 7, 2, 3, 6, 7, 1, 1, 5, 5, 7, 6]

    return mesh3d(x=x, y=y, z=z,
        i=i, j=j, k=k,
        flatshading=true,
        color=color,
        opacity=opc,
        lighting=attr(
            diffuse=0.1,
            specular=1.2,
            roughness=1.0,
        ),
    )
end

"""
    cubes(origin::Vector{<:Real}, side::Real, color::String=""; opc::Real=1)

    Creates a 3D cube mesh centered at the given origin with specified dimensions and color.

    # Arguments
    - `origin::Vector{<:Real}`: A vector of three Reals specifying the center of the cube.
    - `side::Real`: Side length of the cube.
    - `color::String`: A string specifying the color of the cube.
    
    # Keywords
    - `opc`: (optional) A Real specifying the opacity of the cube. Default is 1.                                   
"""
function cubes(origin::Vector{<:Real}, side::Real, color::String=""; opc::Real=1)
    @assert length(origin) == 3
    @assert side > 0

    return cuboids(origin, [side, side, side], color; opc=opc)
end

"""
    squares(origin::Vector{<:Real}, side::Real, mode::String="z", color::String=""; opc::Real=1)
    
    Creates a 2D square mesh centered at the given origin with the specified side length and color.

    # Arguments
    - `origin::Vector{<:Real}`: A vector of three Reals specifying the center of the square.
    - `side::Real`: A Real specifying the side length of the square.
    - `mode`::String: (optional) A string specifying the orientation of the square ("x", "y", or "z"). Default is "z".
    - `color::String`: A string specifying the color of the square.

    
    # Keywords
    - `opc`: (optional) A Real specifying the opacity of the square. Default is 1.
"""
function squares(origin::Vector{<:Real}, side::Real, mode::String="z", color::String=""; opc::Real=1)
    @assert length(origin) == 3

    if color == ""
        @all r g b = round(Int, rand() * 255)
        color = "rgb($r, $g, $b)"
    end

    if mode == "x"
        x = [origin[1], origin[1], origin[1], origin[1]]
        y = [origin[2] - side / 2, origin[2] + side / 2, origin[2] + side / 2, origin[2] - side / 2]
        z = [origin[3] - side / 2, origin[3] - side / 2, origin[3] + side / 2, origin[3] + side / 2]
        i = [0, 2]
        j = [1, 3]
        k = [2, 0]
    elseif mode == "y"
        x = [origin[1] - side / 2, origin[1] - side / 2, origin[1] + side / 2, origin[1] + side / 2]
        y = [origin[2], origin[2], origin[2], origin[2]]
        z = [origin[3] - side / 2, origin[3] + side / 2, origin[3] + side / 2, origin[3] - side / 2]
        i = [0, 2]
        j = [1, 3]
        k = [2, 0]
    else
        x = [origin[1] - side / 2, origin[1] - side / 2, origin[1] + side / 2, origin[1] + side / 2]
        y = [origin[2] - side / 2, origin[2] + side / 2, origin[2] + side / 2, origin[2] - side / 2]
        z = [origin[3], origin[3], origin[3], origin[3]]
        i = [0, 2]
        j = [1, 3]
        k = [2, 0]
    end

    return mesh3d(x=x, y=y, z=z, i=i, j=j, k=k,
        color=color,
        opacity=opc,
        lighting=attr(
            diffuse=0.1,
            specular=1.2,
            roughness=1.0
        ),
    )
end

"""
    ellipsoids(origin::Vector{<:Real}, par::Vector{<:Real}, color::String=""; opc::Real=1, tres=61, pres=31, ah::Real=0)

    Creates a 3D ellipsoid mesh.

    # Arguments
    - `origin::Vector{<:Real}`: The center of the ellipsoid.
    - `par::Vector{<:Real}`: Parameters of the ellipsoid (a, b, c).
    - `color::String`: The color of the ellipsoid.
    
    # Keywords
    - `opc`: The opacity of the ellipsoid. Default is 1.
    - `tres`: The resolution of the mesh grid (theta). Default is 61.
    - `pres`: The resolution of the mesh grid (phi). Default is 31.
    - `ah`: alphahull value. Default is 0.
"""
function ellipsoids(origin::Vector{<:Real}, par::Vector{<:Real}, color::String=""; opc::Real=1, tres=61, pres=31, ah::Real=0)
    @assert length(origin) == 3
    @assert length(par) == 3
    @assert all(par .> 0)
    @assert tres > 0
    @assert pres > 0

    if color == ""
        @all r g b = round(Int, rand() * 255)
        color = "rgb($r, $g, $b)"
    end

    phi = LinRange(0, 360, pres)
    tht = LinRange(0, 180, tres)

    x = sind.(tht) .* cosd.(phi') .* par[1] .+ origin[1]
    y = sind.(tht) .* sind.(phi') .* par[2] .+ origin[2]
    z = cosd.(tht * ones(length(phi))') .* par[3] .+ origin[3]
    x = x[:]
    y = y[:]
    z = z[:]

    return mesh3d(x=x, y=y, z=z,
        alphahull=ah,
        flatshading=true,
        color=color,
        opacity=opc,
        lighting=attr(
            diffuse=0.1,
            specular=2.0,
            roughness=0.5
        ),
    )
end

"""
    spheres(origin::Vector{<:Real}, r::Real, color::String=""; opc::Real=1, tres=60, pres=30, ah::Real=0)

    Creates a 3D sphere mesh.

    # Arguments
    - `origin::Vector{<:Real}`: The center of the ellipsoid.
    - `r::Real`: Radius of the sphere.
    - `color::String`: The color of the ellipsoid.

    # Keywords
    - `opc`: The opacity of the ellipsoid. Default is 1.
    - `tres`: The resolution of the mesh grid (theta). Default is 60.
    - `pres`: The resolution of the mesh grid (phi). Default is 30.
    - `ah`: alphahull value. Default is 0.
"""
function spheres(origin::Vector{<:Real}, r::Real, color::String=""; opc::Real=1, tres=60, pres=30, ah::Real=0)
    @assert length(origin) == 3
    @assert r > 0

    if color == ""
        @all r g b = round(Int, rand() * 255)
        color = "rgb($r, $g, $b)"
    end

    return ellipsoids(origin, [r, r, r], color; opc=opc, tres=tres, pres=pres, ah=ah)
end

"""
    cylinders(origin::Vector{<:Real}, r::Real, h::Real, axis::String="z", color::String=""; opc::Real=1, tres::Int=60, caps::Bool=true)

    Creates a cylinder mesh.

    # Arguments
    - `origin::Vector{<:Real}`: Center point of the cylinder.
    - `r::Real`: Radius.
    - `h::Real`: Height.
    - `axis::String`: Cylinder axis (`"x"`, `"y"`, or `"z"`). Default is `"z"`.
    - `color::String`: Mesh color.

    # Keywords
    - `opc`: Opacity. Default is `1`.
    - `tres`: Circular resolution. Default is `60`.
    - `caps`: Include end caps. Default is `true`.
"""
function cylinders(origin::Vector{<:Real}, r::Real, h::Real, axis::String="z", color::String=""; opc::Real=1, tres::Int=60, caps::Bool=true)
    @assert length(origin) == 3
    @assert r > 0
    @assert h > 0
    @assert tres >= 3
    _validate_axis(axis)

    if color == ""
        color = _random_rgb()
    end

    theta = range(0, 2 * pi, length=tres+1)[1:end-1]
    x = Float64[]
    y = Float64[]
    z = Float64[]

    for th in theta
        push!(x, r * cos(th))
        push!(y, r * sin(th))
        push!(z, -h / 2)
        push!(x, r * cos(th))
        push!(y, r * sin(th))
        push!(z, h / 2)
    end

    i = Int[]
    j = Int[]
    k = Int[]
    for n in 1:tres
        m = n == tres ? 1 : n + 1
        b1 = 2 * n - 1
        t1 = 2 * n
        b2 = 2 * m - 1
        t2 = 2 * m

        push!(i, b1 - 1)
        push!(j, t1 - 1)
        push!(k, t2 - 1)

        push!(i, b1 - 1)
        push!(j, t2 - 1)
        push!(k, b2 - 1)
    end

    if caps
        bcenter = length(x) + 1
        push!(x, 0.0); push!(y, 0.0); push!(z, -h / 2)
        tcenter = length(x) + 1
        push!(x, 0.0); push!(y, 0.0); push!(z, h / 2)

        for n in 1:tres
            m = n == tres ? 1 : n + 1
            b1 = 2 * n - 1
            b2 = 2 * m - 1
            t1 = 2 * n
            t2 = 2 * m

            push!(i, bcenter - 1)
            push!(j, b2 - 1)
            push!(k, b1 - 1)

            push!(i, tcenter - 1)
            push!(j, t1 - 1)
            push!(k, t2 - 1)
        end
    end

    _orient_and_shift!(x, y, z, origin, axis)

    return mesh3d(x=x, y=y, z=z,
        i=i, j=j, k=k,
        color=color,
        opacity=opc,
        lighting=attr(
            diffuse=0.1,
            specular=1.2,
            roughness=1.0,
        ),
    )
end

"""
    cones(origin::Vector{<:Real}, r::Real, h::Real, axis::String="z", color::String=""; opc::Real=1, tres::Int=60, caps::Bool=true)

    Creates a cone mesh.
"""
function cones(origin::Vector{<:Real}, r::Real, h::Real, axis::String="z", color::String=""; opc::Real=1, tres::Int=60, caps::Bool=true)
    return cones(origin, r, 0.0, h, axis, color; opc=opc, tres=tres, caps=caps)
end

"""
    cones(origin::Vector{<:Real}, r1::Real, r2::Real, h::Real, axis::String="z", color::String=""; opc::Real=1, tres::Int=60, caps::Bool=true)

    Creates a cone/frustum mesh. `r1` and `r2` are radii of the bottom and top faces.
"""
function cones(origin::Vector{<:Real}, r1::Real, r2::Real, h::Real, axis::String="z", color::String=""; opc::Real=1, tres::Int=60, caps::Bool=true)
    @assert length(origin) == 3
    @assert r1 >= 0
    @assert r2 >= 0
    @assert r1 + r2 > 0
    @assert h > 0
    @assert tres >= 3
    _validate_axis(axis)

    if color == ""
        color = _random_rgb()
    end

    theta = range(0, 2 * pi, length=tres+1)[1:end-1]
    x = Float64[]
    y = Float64[]
    z = Float64[]
    for th in theta
        push!(x, r1 * cos(th))
        push!(y, r1 * sin(th))
        push!(z, -h / 2)
        push!(x, r2 * cos(th))
        push!(y, r2 * sin(th))
        push!(z, h / 2)
    end

    i = Int[]
    j = Int[]
    k = Int[]
    for n in 1:tres
        m = n == tres ? 1 : n + 1
        b1 = 2 * n - 1
        t1 = 2 * n
        b2 = 2 * m - 1
        t2 = 2 * m

        if r2 > 0
            push!(i, b1 - 1)
            push!(j, t1 - 1)
            push!(k, t2 - 1)
        end
        if r1 > 0
            push!(i, b1 - 1)
            push!(j, t2 - 1)
            push!(k, b2 - 1)
        end
    end

    if caps && r1 > 0
        bcenter = length(x) + 1
        push!(x, 0.0); push!(y, 0.0); push!(z, -h / 2)
        for n in 1:tres
            m = n == tres ? 1 : n + 1
            b1 = 2 * n - 1
            b2 = 2 * m - 1
            push!(i, bcenter - 1)
            push!(j, b2 - 1)
            push!(k, b1 - 1)
        end
    end

    if caps && r2 > 0
        tcenter = length(x) + 1
        push!(x, 0.0); push!(y, 0.0); push!(z, h / 2)
        for n in 1:tres
            m = n == tres ? 1 : n + 1
            t1 = 2 * n
            t2 = 2 * m
            push!(i, tcenter - 1)
            push!(j, t1 - 1)
            push!(k, t2 - 1)
        end
    end

    _orient_and_shift!(x, y, z, origin, axis)

    return mesh3d(x=x, y=y, z=z,
        i=i, j=j, k=k,
        color=color,
        opacity=opc,
        lighting=attr(
            diffuse=0.1,
            specular=1.2,
            roughness=1.0,
        ),
    )
end

"""
    disks(origin::Vector{<:Real}, r::Real, axis::String="z", color::String=""; opc::Real=1, tres::Int=60)

    Creates a circular disk mesh.
"""
function disks(origin::Vector{<:Real}, r::Real, axis::String="z", color::String=""; opc::Real=1, tres::Int=60)
    @assert length(origin) == 3
    @assert r > 0
    @assert tres >= 3
    _validate_axis(axis)

    if color == ""
        color = _random_rgb()
    end

    theta = range(0, 2 * pi, length=tres+1)[1:end-1]
    x = Float64[r * cos(th) for th in theta]
    y = Float64[r * sin(th) for th in theta]
    z = zeros(Float64, tres)
    push!(x, 0.0); push!(y, 0.0); push!(z, 0.0)
    cidx = length(x) - 1

    i = Int[]
    j = Int[]
    k = Int[]
    for n in 1:tres
        m = n == tres ? 1 : n + 1
        push!(i, n - 1)
        push!(j, m - 1)
        push!(k, cidx)
    end

    _orient_and_shift!(x, y, z, origin, axis)

    return mesh3d(x=x, y=y, z=z,
        i=i, j=j, k=k,
        color=color,
        opacity=opc,
        lighting=attr(
            diffuse=0.1,
            specular=1.2,
            roughness=1.0,
        ),
    )
end

"""
    planes(origin::Vector{<:Real}, dims::Vector{<:Real}, axis::String="z", color::String=""; opc::Real=1)

    Creates a rectangular plane mesh where `axis` is the normal direction.
"""
function planes(origin::Vector{<:Real}, dims::Vector{<:Real}, axis::String="z", color::String=""; opc::Real=1)
    @assert length(origin) == 3
    @assert length(dims) == 2
    @assert all(dims .> 0)
    _validate_axis(axis)

    if color == ""
        color = _random_rgb()
    end

    dx = dims[1] / 2
    dy = dims[2] / 2
    x = Float64[-dx, -dx, dx, dx]
    y = Float64[-dy, dy, dy, -dy]
    z = Float64[0, 0, 0, 0]
    i = [0, 2]
    j = [1, 3]
    k = [2, 0]

    _orient_and_shift!(x, y, z, origin, axis)

    return mesh3d(x=x, y=y, z=z,
        i=i, j=j, k=k,
        color=color,
        opacity=opc,
        lighting=attr(
            diffuse=0.1,
            specular=1.2,
            roughness=1.0,
        ),
    )
end

"""
    tori(origin::Vector{<:Real}, R::Real, r::Real, axis::String="z", color::String=""; opc::Real=1, ures::Int=61, vres::Int=31)

    Creates a torus mesh with major radius `R` and minor radius `r`.
"""
function tori(origin::Vector{<:Real}, R::Real, r::Real, axis::String="z", color::String=""; opc::Real=1, ures::Int=61, vres::Int=31)
    @assert length(origin) == 3
    @assert R > 0
    @assert r > 0
    @assert ures >= 3
    @assert vres >= 3
    _validate_axis(axis)

    if color == ""
        color = _random_rgb()
    end

    x = Float64[]
    y = Float64[]
    z = Float64[]

    for iu in 1:ures
        u = 2 * pi * (iu - 1) / ures
        cu = cos(u)
        su = sin(u)
        for iv in 1:vres
            v = 2 * pi * (iv - 1) / vres
            cv = cos(v)
            sv = sin(v)
            push!(x, (R + r * cv) * cu)
            push!(y, (R + r * cv) * su)
            push!(z, r * sv)
        end
    end

    idx(u, v) = (u - 1) * vres + v

    i = Int[]
    j = Int[]
    k = Int[]
    for iu in 1:ures
        nu = iu == ures ? 1 : iu + 1
        for iv in 1:vres
            nv = iv == vres ? 1 : iv + 1
            a = idx(iu, iv)
            b = idx(nu, iv)
            c = idx(nu, nv)
            d = idx(iu, nv)

            push!(i, a - 1); push!(j, b - 1); push!(k, c - 1)
            push!(i, a - 1); push!(j, c - 1); push!(k, d - 1)
        end
    end

    _orient_and_shift!(x, y, z, origin, axis)

    return mesh3d(x=x, y=y, z=z,
        i=i, j=j, k=k,
        color=color,
        opacity=opc,
        lighting=attr(
            diffuse=0.1,
            specular=1.2,
            roughness=1.0,
        ),
    )
end

"""
    lines(pt1::Vector{<:Real}, pt2::Vector{<:Real}, color::String; opc::Real=1, style="")

    Creates a 3D line between two points.

    # Arguments
    - `pt1::Vector{<:Real}`: Starting point of the line.
    - `pt2::Vector{<:Real}`: Ending point of the line.
    - `color::String`: The color of the line.

    # Keywords
    - `opc`: The opacity of the line. Default is 1.
    - `style`: The line style (e.g., "solid", "dash"). Default is "".
"""
function lines(pt1::Vector{<:Real}, pt2::Vector{<:Real}, color::String=""; opc::Real=1, style="")
    @assert length(pt1) == 3
    @assert length(pt2) == 3

    if color == ""
        @all r g b = round(Int, rand() * 255)
        color = "rgb($r, $g, $b)"
    end

    x = [pt1[1], pt2[1]]
    y = [pt1[2], pt2[2]]
    z = [pt1[3], pt2[3]]

    return scatter3d(x=x, y=y, z=z,
        mode="lines",
        line=attr(
            color=color,
            dash=style,
        ),
        showlegend=false,
        flatshading=true,
        opacity=opc,
    )
end

"""
    polygons(pts::Vector, color::String; opc::Real=1, ah::Real=0)

    Creates a polygon mesh from a set of points (form around the mid point of the set of points).

    # Arguments
    - `pts::Vector`: List of points defining the polygon.
    - `color::String`: The color of the polygon.

    # Keywords
    - `opc`: The opacity of the polygon. Default is 1.
    - `ah`: alphahull value. Default is 0.
"""
function polygons(pts::Vector, color::String=""; opc::Real=1, ah::Real=0)
    _validate_points3d(pts)
    @assert length(pts) >= 3

    if color == ""
        color = _random_rgb()
    end

    x = Float64[]
    y = Float64[]
    z = Float64[]
    i = Int[]
    j = Int[]
    k = Int[]
    _append_polygon_mesh!(x, y, z, i, j, k, pts, 0)

    return mesh3d(x=x, y=y, z=z,
        i=i, j=j, k=k,
        alphahull=ah,
        color=color,
        opacity=opc,
        lighting=attr(
            diffuse=0.1,
            specular=1.2,
            roughness=1.0
        ),
    )
end

"""
    polygons(pts::Vector, ng::Int, color::String; opc::Real=1, ah::Real=0)

    Creates a group of polygons from a set of points and a specified number of vertices per polygon.

    # Arguments
    - `pts::Vector`: List of points defining the mesh.
    - `ng::Int`: Number of vertices per polygon.
    - `color::String`: The color of the mesh.

    # Keywords
    - `opc`: The opacity of the mesh. Default is 1.
    - `ah`: alphahull value.
"""
function polygons(pts::Vector, ng::Int, color::String=""; opc::Real=1, ah::Real=0)
    _validate_points3d(pts)
    @assert ng > 0
    @assert length(pts) >= ng
    @assert mod(length(pts), ng) == 0 "Point count must be a multiple of ng."

    if color == ""
        color = _random_rgb()
    end

    x = Float64[]
    y = Float64[]
    z = Float64[]
    i = Int[]
    j = Int[]
    k = Int[]
    offset = 0

    N = length(pts)
    Ng = floor(Int, N / ng)
    for p = 1:Ng
        ptsg = []
        for m = 1:ng
            push!(ptsg, pts[(p-1)*ng+m])
        end
        offset = _append_polygon_mesh!(x, y, z, i, j, k, ptsg, offset)
    end

    return mesh3d(x=x, y=y, z=z,
        i=i, j=j, k=k,
        alphahull=ah,
        color=color,
        opacity=opc,
        lighting=attr(
            diffuse=0.1,
            specular=1.2,
            roughness=1.0
        ),
    )
end

"""
    stlmesh(filepath::String, color::String=""; opc::Real=1)

    Loads an STL file (binary or ASCII) and returns a mesh3d trace.

    # Arguments
    - `filepath::String`: Path to the STL file.
    - `color::String`: The color of the mesh.

    # Keywords
    - `opc`: The opacity of the mesh. Default is 1.
"""
function stlmesh(filepath::String, color::String=""; opc::Real=1)
    @assert isfile(filepath) "File not found: $filepath"

    if color == ""
        color = _random_rgb()
    end

    x, y, z, i, j, k = open(filepath) do io
        if _is_stl_binary(filepath)
            _parse_stl_binary(io)
        else
            _parse_stl_ascii(io)
        end
    end
    @assert length(x) > 0 "No triangles found in STL file."

    return mesh3d(x=x, y=y, z=z,
        i=i, j=j, k=k,
        flatshading=true,
        color=color,
        opacity=opc,
        lighting=attr(
            diffuse=0.1,
            specular=1.2,
            roughness=1.0
        ),
    )
end

"""
    gtrans!(geo::GenericTrace, dis::Vector{<:Real})

    Translates a 3D geometry by a specified displacement vector.

    # Arguments
    - `geo::GenericTrace`: The geometry to translate.
    - `dis::Vector{<:Real}`: A vector of three Reals specifying the translation distances for the x, y, and z axes.
"""
function gtrans!(geo::GenericTrace, dis::Vector{<:Real})
    @assert length(dis) == 3
    @inbounds for n in eachindex(geo.x)
        geo.x[n] += dis[1]
        geo.y[n] += dis[2]
        geo.z[n] += dis[3]
    end
    return nothing
end

"""
    gtrans!(geos::AbstractVector{<:GenericTrace}, dis::Vector{<:Real})

    Translates each geometry in a vector by a specified displacement vector.
"""
function gtrans!(geos::AbstractVector{<:GenericTrace}, dis::Vector{<:Real})
    @assert length(dis) == 3
    for geo in geos
        gtrans!(geo, dis)
    end
    return nothing
end

"""
    grot!(geo::GenericTrace, rotang::Vector{<:Real}, center::Vector{<:Real}=[0])

    Rotates a 3D geometry around a specified center point. (Tait–Bryan rotation)

    # Arguments
    - `geo::GenericTrace`: The 3D geometry to be rotated, which must have `x`, `y`, and `z` coordinates.
    - `rotang::Vector{<:Real}`: A vector of three Tait–Bryan rotation angles in degrees  for rotations around the x, y, and z axes respectively.
    - `center::Vector{<:Real}`: The center point of rotation. Default is `[0]`, which means the rotation center will be set at the geometric center of the object.
"""
function grot!(geo::GenericTrace, rotang::Vector{<:Real}, center::Vector{<:Real}=[0])
    @assert length(rotang) == 3

    alpha = rotang[1]
    beta = rotang[2]
    gama = rotang[3]

    Rx = [1 0 0;
        0 cosd(alpha) -sind(alpha);
        0 sind(alpha) cosd(alpha)]
    Ry = [cosd(beta) 0 sind(beta);
        0 1 0;
        -sind(beta) 0 cosd(beta)]
    Rz = [cosd(gama) -sind(gama) 0;
        sind(gama) cosd(gama) 0;
        0 0 1]

    R = Rz * Ry * Rx

    pos = []
    for n in eachindex(geo.x)
        push!(pos, [geo.x[n], geo.y[n], geo.z[n]])
    end

    if center == [0] # rotation center set at the geometry center
        center = sum(pos) ./ length(pos)
    else
        @assert length(center) == 3
    end

    @inbounds for n in eachindex(geo.x)
        vec = R * (pos[n] - center)
        geo.x[n] = vec[1] + center[1]
        geo.y[n] = vec[2] + center[2]
        geo.z[n] = vec[3] + center[3]
    end
    return nothing
end

"""
    grot!(geos::AbstractVector{<:GenericTrace}, rotang::Vector{<:Real}, center::Vector{<:Real}=[0])

    Rotates each geometry in a vector using Tait-Bryan angles.
"""
function grot!(geos::AbstractVector{<:GenericTrace}, rotang::Vector{<:Real}, center::Vector{<:Real}=[0])
    for geo in geos
        grot!(geo, rotang, center)
    end
    return nothing
end

"""
    grot!(geo::GenericTrace, ang::Real, axis::Vector{<:Real}, origin::Vector{<:Real}=[0])

    Rotates the geometry by the specified angle `ang` around the axis `axis` and origin `origin`.
    
    # Arguments
    - `geo::GenericTrace`: The geometry to be rotated.
    - `ang::Real`: The rotation angle.
    - `axis::Vector{<:Real}`: The rotation axis.
    - `origin::Vector{<:Real}=[0]`: The rotation origin. Defaults to the center of the geometry if not specified.
"""
function grot!(geo::GenericTrace, ang::Real, axis::Vector{<:Real}, origin::Vector{<:Real}=[0])
    @assert length(axis) == 3
    @assert norm(axis) > 0

    pos = []
    for n in eachindex(geo.x)
        push!(pos, [geo.x[n], geo.y[n], geo.z[n]])
    end

    axis_norm = axis ./ norm(axis)
    vrot = similar(pos)

    if origin == [0] # rotation center set at the geometry center
        origin = sum(pos) ./ length(pos)
    else
        @assert length(origin) == 3
    end

    v = similar(origin)
    for n in eachindex(vrot)
        v .= (pos[n] .- origin)
        vrot[n] = cosd(ang) * v + sind(ang) * cross(axis_norm, v) + (1 - cosd(ang)) * dot(axis_norm, v) * axis_norm
        pos[n] = vrot[n] .+ origin
    end

    geo.x = getindex.(pos, 1)
    geo.y = getindex.(pos, 2)
    geo.z = getindex.(pos, 3)
    return nothing
end

"""
    grot!(geos::AbstractVector{<:GenericTrace}, ang::Real, axis::Vector{<:Real}, origin::Vector{<:Real}=[0])

    Rotates each geometry in a vector around a given axis and origin.
"""
function grot!(geos::AbstractVector{<:GenericTrace}, ang::Real, axis::Vector{<:Real}, origin::Vector{<:Real}=[0])
    @assert norm(axis) > 0
    for geo in geos
        grot!(geo, ang, axis, origin)
    end
    return nothing
end

"""
    sort_pts(pts::Vector)

    Sorts points in place based on their angular position relative to the centroid.

    # Arguments
    - `pts::Vector`: List of points to be sorted.
"""
function sort_pts(pts::Vector)
    @assert all(length.(pts) .== 3)

    N = length(pts)
    ang = zeros(length(pts))
    mid = zeros(3)
    for n in eachindex(pts)
        for m in 1:3
            mid[m] += pts[n][m]
        end
    end
    mid .= mid ./ N

    c = collect(combinations(1:N, 3))
    thtr = 0
    phir = 0

    vec = similar(mid)
    for n in eachindex(c)
        vec .= cross(pts[c[n][2]] .- pts[c[n][1]], pts[c[n][3]] .- pts[c[n][1]])
        if norm(vec) == 0
            continue
        else
            vec .= vec ./ norm(vec)
            thtr = acosd(vec[3])
            phir = atand(vec[2], vec[1])
            break
        end
    end

    Ry = [
        cosd(thtr) 0 -sind(thtr);
        0 1 0;
        sind(thtr) 0 cosd(thtr)
    ]
    Rz = [
        cosd(phir) sind(phir) 0;
        -sind(phir) cosd(phir) 0;
        0 0 1
    ]
    R = Ry * Rz
    pts_rot = similar(pts)
    mid_rot = R * mid
    for n in eachindex(pts)
        pts_rot[n] = R * pts[n]
        ang[n] = atan(pts_rot[n][2] - mid_rot[2], pts_rot[n][1] - mid_rot[1])
    end
    pts_rot .= pts[sortperm(ang)]

    return pts_rot
end

"""
    sort_pts!(pts::Vector)

    Sorts points in place based on their angular position relative to the centroid.

    # Arguments
    - `pts::Vector`: List of points to be sorted.
"""
function sort_pts!(pts::Vector)
    @assert all(length.(pts) .== 3)

    N = length(pts)
    ang = zeros(length(pts))
    mid = zeros(3)
    for n in eachindex(pts)
        for m in 1:3
            mid[m] += pts[n][m]
        end
    end
    mid .= mid ./ N

    c = collect(combinations(1:N, 3))
    thtr = 0
    phir = 0
    vec = similar(mid)
    for n in eachindex(c)
        vec .= cross(pts[c[n][2]] .- pts[c[n][1]], pts[c[n][3]] .- pts[c[n][1]])
        if norm(vec) == 0
            continue
        else
            vec .= vec ./ norm(vec)
            thtr = acosd(vec[3])
            phir = atand(vec[2], vec[1])
            break
        end
    end

    Ry = [
        cosd(thtr) 0 -sind(thtr);
        0 1 0;
        sind(thtr) 0 cosd(thtr)
    ]
    Rz = [
        cosd(phir) sind(phir) 0;
        -sind(phir) cosd(phir) 0;
        0 0 1
    ]
    R = Ry * Rz
    pts_rot = similar(pts)
    mid_rot = R * mid
    for n in eachindex(pts)
        pts_rot[n] = R * pts[n]
        ang[n] = atan(pts_rot[n][2] - mid_rot[2], pts_rot[n][1] - mid_rot[1])
    end
    pts .= pts[sortperm(ang)]
    return nothing
end

"""
    add_ref_axes!(plt::Union{Plot, SyncPlot}, origin::Vector{<:Real}=[0, 0, 0], r::Real=1)

    Adds reference axes (x, y, z) to a plot.

    # Arguments
    - `plt::Union{Plot, SyncPlot}`: The plot to which the axes will be added.
    - `origin::Vector{<:Real}`: The origin point of the axes.
    - `r::Real`: The length of the reference axes.
"""
function add_ref_axes!(plt::PlotContainer, origin::Vector{<:Real}=[0, 0, 0], r::Real=1)
    @assert length(origin) == 3
    @assert r > 0

    cx = cone(x=[r + origin[1]], y=[origin[2]], z=[origin[3]], u=[r / 10], v=[0], w=[0],
        colorscale=[[0, "red"], [1, "red"]],
        showscale=false)
    cy = cone(x=[origin[1]], y=[r + origin[2]], z=[origin[3]], u=[0], v=[r / 10], w=[0],
        colorscale=[[0, "green"], [1, "green"]],
        showscale=false)
    cz = cone(x=[origin[1]], y=[origin[2]], z=[r + origin[3]], u=[0], v=[0], w=[r / 10],
        colorscale=[[0, "blue"], [1, "blue"]],
        showscale=false)
    lx = scatter3d(x=[origin[1], r + origin[1]], y=[origin[2], origin[2]], z=[origin[3], origin[3]],
        line=attr(color="red"),
        mode="lines",
        showlegend=false)
    ly = scatter3d(x=[origin[1], origin[1]], y=[origin[2], r + origin[2]], z=[origin[3], origin[3]],
        line=attr(color="green"),
        mode="lines",
        showlegend=false)
    lz = scatter3d(x=[origin[1], origin[1]], y=[origin[2], origin[2]], z=[origin[3], r + origin[3]],
        line=attr(color="blue"),
        mode="lines",
        showlegend=false)
    addtraces!(plt, cx)
    addtraces!(plt, cy)
    addtraces!(plt, cz)
    addtraces!(plt, lx)
    addtraces!(plt, ly)
    addtraces!(plt, lz)
    relayout!(plt, scene=attr(
        annotations=[
            attr(
                showarrow=false,
                x=origin[1] + 1.1 * r,
                y=origin[2],
                z=origin[3],
                text="x",
                font=attr(color="red")
            ),
            attr(
                showarrow=false,
                x=origin[1],
                y=origin[2] + 1.1 * r,
                z=origin[3],
                text="y",
                font=attr(color="green")
            ),
            attr(
                showarrow=false,
                x=origin[1],
                y=origin[2],
                z=origin[3] + 1.1 * r,
                text="z",
                font=attr(color="blue")
            ),
        ]
    ))
    return nothing
end

"""
    add_ref_axes!(plt::Union{Plot, SyncPlot}, origin::Vector{<:Real}, r::Vector{<:Real})

    Adds reference axes (x, y, z) to a plot.

    # Arguments
    - `plt::Union{Plot, SyncPlot}`: The plot to which the axes will be added.
    - `origin::Vector{<:Real}`: The origin point of the axes.
    - `r::Vector{<:Real}`: The lengths of the reference axes.
"""
function add_ref_axes!(plt::PlotContainer, origin::Vector{<:Real}, r::Vector{<:Real})
    @assert length(origin) == 3
    @assert length(r) == 3
    @assert all(r .> 0)

    csize = minimum([r[1] / 10, r[2] / 10, r[3] / 10])
    cx = cone(x=[r[1] + origin[1]], y=[origin[2]], z=[origin[3]], u=[csize], v=[0], w=[0],
        colorscale=[[0, "red"], [1, "red"]],
        showscale=false)
    cy = cone(x=[origin[1]], y=[r[2] + origin[2]], z=[origin[3]], u=[0], v=[csize], w=[0],
        colorscale=[[0, "green"], [1, "green"]],
        showscale=false)
    cz = cone(x=[origin[1]], y=[origin[2]], z=[r[3] + origin[3]], u=[0], v=[0], w=[csize],
        colorscale=[[0, "blue"], [1, "blue"]],
        showscale=false)
    lx = scatter3d(x=[origin[1], r[1] + origin[1]], y=[origin[2], origin[2]], z=[origin[3], origin[3]],
        line=attr(color="red"),
        mode="lines",
        showlegend=false)
    ly = scatter3d(x=[origin[1], origin[1]], y=[origin[2], r[2] + origin[2]], z=[origin[3], origin[3]],
        line=attr(color="green"),
        mode="lines",
        showlegend=false)
    lz = scatter3d(x=[origin[1], origin[1]], y=[origin[2], origin[2]], z=[origin[3], r[3] + origin[3]],
        line=attr(color="blue"),
        mode="lines",
        showlegend=false)
    addtraces!(plt, cx)
    addtraces!(plt, cy)
    addtraces!(plt, cz)
    addtraces!(plt, lx)
    addtraces!(plt, ly)
    addtraces!(plt, lz)
    relayout!(plt, scene=attr(
        annotations=[
            attr(
                showarrow=false,
                x=origin[1] + 1.1 * r[1],
                y=origin[2],
                z=origin[3],
                text="x",
                font=attr(color="red")
            ),
            attr(
                showarrow=false,
                x=origin[1],
                y=origin[2] + 1.1 * r[2],
                z=origin[3],
                text="y",
                font=attr(color="green")
            ),
            attr(
                showarrow=false,
                x=origin[1],
                y=origin[2],
                z=origin[3] + 1.1 * r[3],
                text="z",
                font=attr(color="blue")
            ),
        ]
    ))
    return nothing
end

"""
    add_arrows!(plt::Union{Plot, SyncPlot}, origin::Vector{<:Real}, dir::Vector{<:Real}, len::Real=1.0, color::String=""; opc::Real=1, endpoint::Bool=true, asize::Real=len)
    
    Creates a 3D arrow starting from a point and pointing in a given direction.

    # Arguments
    - `plt::Union{Plot, SyncPlot}`: The plot to which the axes will be added.
    - `origin::Vector{<:Real}`: The starting point of the arrow.
    - `dir::Vector{<:Real}`: The direction vector of the arrow.
    - `len::Real`: length of the arrow
    - `color::String`: The color of the arrow.

    # Keywords
    - `opc`: The opacity of the arrow. Default is 1.
    - `asize` Size factor of the arrow cone. Default is `len`.
"""
function add_arrows!(plt::PlotContainer, origin::Vector{<:Real}, dir::Vector{<:Real}, len::Real=1.0, color::String=""; opc::Real=1, endpoint::Bool=true, asize::Real=len)
    @assert length(origin) == 3
    @assert length(dir) == 3
    @assert len > 0
    @assert asize > 0
    @assert norm(dir) > 0 "Direction vector must be non-zero."

    if color == ""
        color = _random_rgb()
    end

    dir_norm = dir / norm(dir) * len
    lwidth = max(1.0, 4 * asize / len)

    if endpoint == true
        c = cone(x=[origin[1] + dir_norm[1]], y=[origin[2] + dir_norm[2]], z=[origin[3] + dir_norm[3]], u=[dir_norm[1]/len*asize], v=[dir_norm[2]]/len*asize, w=[dir_norm[3]/len*asize],
            colorscale=[[0, color], [1, color]],
            opacity=opc,
            showscale=false)
        l = scatter3d(x=[origin[1], origin[1] + dir_norm[1]], y=[origin[2], origin[2] + dir_norm[2]], z=[origin[3], origin[3] + dir_norm[3]],
            line=attr(color=color, width=lwidth),
            mode="lines",
            opacity=opc,
            showlegend=false)
    else # origin at center
        c = cone(x=[origin[1] + dir_norm[1] / 2], y=[origin[2] + dir_norm[2] / 2], z=[origin[3] + dir_norm[3] / 2], u=[dir_norm[1]/len*asize], v=[dir_norm[2]/len*asize], w=[dir_norm[3]/len*asize],
            colorscale=[[0, color], [1, color]],
            opacity=opc,
            showscale=false)
        l = scatter3d(x=[origin[1] - dir_norm[1] / 2, origin[1] + dir_norm[1] / 2], y=[origin[2] - dir_norm[2] / 2, origin[2] + dir_norm[2] / 2], z=[origin[3] - dir_norm[3] / 2, origin[3] + dir_norm[3] / 2],
            line=attr(color=color, width=lwidth),
            mode="lines",
            opacity=opc,
            showlegend=false)
    end

    addtraces!(plt, c)
    addtraces!(plt, l)

    return nothing
end

"""
    add_text!(plt::Union{Plot, SyncPlot}, origin::Vector{<:Real}, text::String, color::String="")

    Add text to plot. 

    # Arguments
    - plt::Union{Plot, SyncPlot}: Plot to add text.
    - origin::Vector{<:Real}: origin of the text.
    - text::String:: text to be added.
    - color::String=: color of the text.
"""
function add_text!(plt::PlotContainer, origin::Vector{<:Real}, text::String, color::String="")

    if color == ""
        @all r g b = round(Int, rand() * 255)
        color = "rgb($r, $g, $b)"
    end

    addtraces!(plt, scatter3d(x=[origin[1]], y=[origin[2]], z=[origin[3]],
        mode="text", text=[text], textposition="middle center",
        textfont=attr(
            color=color,
        ),
        showlegend=false,
    ))
    return nothing
end

"""
    blank_layout()

    Return blank layout.
"""
function blank_layout()
    layout = Layout(
        scene=attr(
            xaxis=attr(
                visible=false,
                showgrid=false
            ),
            yaxis=attr(
                visible=false,
                showgrid=false
            ),
            zaxis=attr(
                visible=false,
                showgrid=false
            ),
            aspectmode="data",
        ),
    )
    return layout
end

"""
    set_view!(plt::Union{Plot, SyncPlot}, az::Real, el::Real, r::Real=1)

    Set az/el (deg) view of the plot.

    # Arguments
    - plt::Union{Plot, SyncPlot}: plot to be modified.
    - az::Real: az value.
    - el::Real: el value.
    - r::Real=1: distance.
"""
function set_view!(plt::PlotContainer, az::Real, el::Real, r::Real=1.25 * sqrt(3))
    if el == 90
        el -= 0.01
    elseif el == -90
        el += 0.01
    end
    x = r * cosd(el) * cosd(az)
    y = r * cosd(el) * sind(az)
    z = r * sind(el)
    relayout!(plt, scene_camera=attr(eye=attr(x=x, y=y, z=z)))
    return nothing
end
