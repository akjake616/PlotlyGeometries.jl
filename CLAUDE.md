# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

PlotlyGeometries.jl is a Julia package for constructing and manipulating 3D geometry traces for Plotly visualization. It builds on PlotlySupply to provide geometry constructors (cuboids, spheres, cylinders, polygons, etc.), in-place transforms (translation, rotation), and plot helpers (axes, arrows, text, camera).

## Development Commands

```bash
# Run tests
julia --project -e 'using Pkg; Pkg.test()'

# Start a REPL with the package loaded
julia --project -e 'using PlotlyGeometries'

# Build/instantiate dependencies
julia --project -e 'using Pkg; Pkg.instantiate()'
```

## Architecture

The package is essentially two files:
- **src/PlotlyGeometries.jl** — Module definition, imports, and exports
- **src/api.jl** — All implementation (~1400 lines), organized into sections:
  1. **Internal helpers** (top) — `_GEOM_EPS`, `_orient_point`, `_orient_and_shift!`, 2D projection for triangulation, ear-clipping triangulation algorithm
  2. **Geometry constructors** — Each returns a PlotlySupply trace (`mesh3d` or `scatter3d`): `cuboids`, `cubes`, `squares`, `ellipsoids`, `spheres`, `cylinders`, `cones`, `frustums`, `disks`, `planes`, `tori`, `polygons`, `lines`
  3. **Transforms** — In-place mutation via `gtrans!` (translation) and `grot!` (rotation); both accept single traces or vectors of traces
  4. **Plot helpers** — `add_ref_axes!`, `add_arrows!`, `add_text!`, `blank_layout`, `set_view!`

## Key Conventions

- All angles are in **degrees** (not radians)
- Rotation uses Tait-Bryan order (Z-Y-X) or axis-angle via Rodrigues' formula
- Coordinates are `[x, y, z]` vectors
- Constructors accept an axis orientation string (`"x"`, `"y"`, or `"z"`) for shapes like cylinders and cones
- Color is random RGB if not specified; opacity defaults to 1.0
- Input validation uses `@assert` (throws `AssertionError`)
- `PlotContainer = Union{Plot, SyncPlot}` is the accepted figure type for plot-mutating functions
- Resolution keywords: `tres` (circumferential), `ures`/`vres` (toroidal)

## Dependencies

- **PlotlySupply** (>=1.6) — Plotly backend providing `Plot`, `SyncPlot`, `mesh3d`, `scatter3d`, etc.
- **BatchAssign** — `@all` macro for batch field assignment
- **Combinatorics** — `combinations()` used in polygon point sorting
- **LinearAlgebra** — Cross products, norms, matrix operations

## Testing

Tests are in `test/runtests.jl` with four `@testset` groups: base geometry + plot mutation, new primitives, concave polygon triangulation, and group transforms + edge guards. Tests validate trace types, coordinate values (atol=1e-9), plot data structure, and assertion errors.
