# Finite Element Simulation of Incompressible Stokes Flow (P2-P1 Taylor-Hood Elements)

**Author:** Louen MARX

## Description

This project implements a **Finite Element Method (FEM)** solver in MATLAB for simulating **steady-state, viscous, incompressible fluid flows** governed by the **Stokes equations**. The Stokes equations model fluid behavior at low Reynolds numbers by neglecting convective terms in favor of diffusion:

```
-ν Δu + ∇p = f    in Ω
    div(u) = 0     in Ω
```

where:
- `u = (u₁, u₂)` is the velocity field (2D)
- `p` is the pressure
- `ν` is the kinematic viscosity
- `f` is the external body force

The project covers the full FEM pipeline: mesh reading (Gmsh format), element matrix assembly, global system construction with block saddle-point structure, boundary condition enforcement via pseudo-elimination, solving, error analysis, and visualization.

---

## Numerical Methodology

### 1. Taylor-Hood (P2-P1) Elements

To ensure numerical stability and satisfy the **inf-sup (Brezzi-Babuška) condition**, the project uses the **Taylor-Hood** element pair:

| Field | Element | Nodes per triangle | Description |
|-------|---------|-------------------|-------------|
| Velocity (u₁, u₂) | P2 (quadratic) | 6 (3 vertices + 3 edge midpoints) | Higher-order approximation for velocity |
| Pressure (p) | P1 (linear) | 3 (vertices only) | Standard linear interpolation |

The P2-P1 pairing is the standard stable choice for Stokes problems. The project also includes a **P1-P1 implementation** for comparison (which is known to be unstable and used as a pedagogical counterexample).

### 2. Block Saddle-Point System

The coupled velocity-pressure system is assembled into a global **saddle-point matrix** structured by blocks:

```
| νK   0    E |   | u₁ |   | 0 |
| 0    νK   F | × | u₂ | = | 0 |
| Eᵀ   Fᵀ   0 |   | p  |   | 0 |
```

Where:
- `K` = stiffness matrix (viscous diffusion, P2)
- `M` = mass matrix (for error computation, P2)
- `E, F` = velocity-pressure coupling blocks (gradient/divergence terms)
- The `(3,3)` zero block reflects the incompressibility constraint

### 3. Quadrature

Element matrices are computed using **Gauss-Legendre quadrature** on the reference triangle:
- **6-point quadrature** (degree 4) for mass matrices (P2)
- **3-point quadrature** (degree 2) for stiffness matrices (P2)
- Affine transformations map from reference to physical elements

### 4. Boundary Conditions

Boundary conditions are enforced via **pseudo-elimination** (also called penalty-free Dirichlet):

- **Dirichlet (velocity):** Prescribed velocity at walls and inlet. Rows and columns are zeroed out, diagonal set to 1, and RHS set to the boundary value. The RHS for interior nodes is corrected by subtracting the contribution of known boundary values.
- **Neumann (natural):** Applied at the outlet (x = x_max). No explicit enforcement needed — natural boundary conditions are automatically satisfied by the variational formulation.
- **Pressure gauge:** Since pressure is defined up to a constant in the Stokes equations, `p = 0` is fixed at one point to ensure uniqueness.

### 5. Convergence Analysis

The code performs **mesh convergence studies** across multiple mesh sizes (h = 0.2, 0.1, 0.05, 0.025). Errors are computed in:
- **L² norm** (relative): measures the overall error in the solution
- **H¹ semi-norm** (relative): measures the error in the gradient (derivatives)

Log-log regression is used to extract the **observed order of convergence** and compared against theoretical predictions for P2 elements (O(h³) in L², O(h²) in H¹).

---

## Test Cases

### Poiseuille Flow (Validation)

The main validation case is the **Poiseuille flow** in a rectangular channel `[0,3] × [0,2]`:

```
u₁(x,y) = (2 - y) · y        (parabolic profile)
u₂(x,y) = 0
p(x,y)  = -2(x - 3)
```

Boundary conditions:
- **Inlet (x = 0):** Poiseuille profile `u₁ = (2-y)y, u₂ = 0`
- **Walls (y = 0, y = 2):** No-slip `u = 0`
- **Outlet (x = 3):** Natural Neumann (free outflow)

The flow rate `Q` is computed numerically using **trapezoidal quadrature** across the inlet section and compared to the exact value `Q = 4/3`.

### Channel with Cavity (Backward-Facing Step)

A more complex geometry with a cavity attached to the bottom of the channel:

```
Domain:
  - Main channel: x ∈ [0, 3], y ∈ [0, 1]
  - Cavity:       x ∈ [0, 1], y ∈ [-0.5, 0]
```

Boundary conditions:
- **Inlet (left, y ≥ 0):** Poiseuille profile `u₁ = (1-y)y, u₂ = 0`
- **Cavity + walls:** No-slip `u = 0`
- **Outlet (right):** Natural Neumann
- **Pressure:** Fixed `p = 0` at one point

This test case produces a **recirculation zone** inside the cavity, visible in the velocity vector field.

---

## Project Structure

### Main Scripts

| File | Description |
|------|-------------|
| `principal_stokes.m` | Main Stokes solver — Poiseuille flow validation with convergence analysis across 4 meshes |
| `principal_stokes1.m` | Stokes solver for the **channel with cavity** geometry |
| `principalstokesP1.m` | Stokes solver using **P1-P1 elements** (unstable, for comparison) |
| `principal_dirichlet_p2.m` | Scalar P2 FEM for `-Δu + u = f` with **Dirichlet BC** + convergence study |
| `principal_neumann_p2.m` | Scalar P2 FEM for `-Δu + u = f` with **Neumann BC** + convergence study + 1D cross-sections |

### Element Matrix Functions

| File | Description |
|------|-------------|
| `matM_elem_p2.m` | P2 mass matrix on a triangle (6×6) |
| `matK_elem_p2.m` | P2 stiffness matrix on a triangle (6×6) |
| `matE_elem.m` | Velocity-pressure coupling block `∂φ_j/∂x · ψ_i` (P2-P1, 6×3) |
| `matF_elem.m` | Velocity-pressure coupling block `∂φ_j/∂y · ψ_i` (P2-P1, 6×3) |
| `matM_elemP1_P1.m` | P1 mass matrix (3×3) |
| `matk_elemP1_P1.m` | P1 stiffness matrix (3×3) |
| `matE_elemP1_P1.m` | P1 velocity-pressure coupling block `∂φ_j/∂x · ψ_i` (3×3) |
| `matF_elemP1_P1.m` | P1 velocity-pressure coupling block `∂φ_j/∂y · ψ_i` (3×3) |

### Boundary Condition Functions

| File | Description |
|------|-------------|
| `elimine_stokes.m` | Pseudo-elimination for Stokes — Poiseuille flow geometry (rectangular channel) |
| `elimine_stokes1.m` | Pseudo-elimination for Stokes — channel with cavity geometry |
| `elimine.m` | Pseudo-elimination for scalar Dirichlet problem (P2) |
| `elimine2.m` | Pseudo-elimination for scalar homogeneous Dirichlet problem |

### Source Term / Boundary Functions

| File | Description |
|------|-------------|
| `f.m` | Right-hand side `f(x,y)` for Neumann problem |
| `f1.m` | Right-hand side `f(x,y)` for Dirichlet problem |
| `g.m` | Dirichlet boundary data `g(x,y)` |
| `g1.m` | Velocity BC `g₁(x,y)` (x-component) for Stokes |
| `g2.m` | Velocity BC `g₂(x,y)` (y-component) for Stokes |
| `g11.m` | Extended velocity BC with full Poiseuille profile |

### Visualization Functions

| File | Description |
|------|-------------|
| `affiche.m` | Plot scalar field on P1 triangular mesh |
| `affiche_ordre1.m` | Plot scalar field with P1 (linear) elements |
| `affiche_ordre2.m` | Plot scalar field with P2 (quadratic) elements |
| `affichemaillage.m` | Display the P1 mesh structure |
| `affichemaillage_ordre2.m` | Display the P2 mesh structure |

### Mesh I/O

| File | Description |
|------|-------------|
| `lecture_msh.m` | Read Gmsh `.msh` files (order 1, P1 nodes) |
| `lecture_msh_ordre2.m` | Read Gmsh `.msh` files (order 2, P2 nodes with edge midpoints) |

### Mesh Files (Gmsh)

| File | Mesh Size | Description |
|------|-----------|-------------|
| `geomRectangle.msh` | h = 0.2 | Coarse rectangular mesh |
| `geomRectangle2.msh` | h = 0.1 | Medium rectangular mesh |
| `geomRectangle3.msh` | h = 0.05 | Fine rectangular mesh |
| `geomRectangle4.msh` | h = 0.025 | Very fine rectangular mesh |
| `geomCarreP1_P1.msh` | h = 0.1 | Square mesh for P1-P1 test |
| `geomRectangle_partie3.msh` | h = 0.2 | Stokes test mesh (coarse) |
| `geomRectangle_partie3-2.msh` | h = 0.1 | Stokes test mesh (medium) |
| `geomRectangle_partie3-3.msh` | h = 0.05 | Stokes test mesh (fine) |
| `geomRectangle_partie3-4.msh` | h = 0.0375 | Stokes test mesh (very fine) |
| `geomRectangle_partie3-5.msh` | — | Channel with cavity mesh |
| `*.geo` files | — | Gmsh geometry scripts (for mesh generation) |

---

## How to Run

### Requirements

- **MATLAB** (R2020 or later recommended)
- No external toolboxes required (only built-in MATLAB functions are used)

### Execution

Open MATLAB, navigate to the project directory, and run:

```matlab
% Scalar P2 problems:
principal_dirichlet_p2    % Dirichlet problem with convergence analysis
principal_neumann_p2      % Neumann problem with convergence analysis

% Stokes flow:
principal_stokes          % Poiseuille flow (P2-P1 Taylor-Hood)
principal_stokes1         % Channel with cavity (P2-P1 Taylor-Hood)
principalstokesP1         % P1-P1 Stokes (unstable comparison)
```

### Output

Each script produces:
- **Console output:** Error norms (L², H¹), convergence rates, flow rates
- **MATLAB figures:** Velocity fields, pressure maps, error distributions, quiver plots, convergence curves (log-log)

---

## Expected Convergence Results

### P2 Elements (Dirichlet/Neumann scalar problems)

| Norm | Theoretical Order | Expected Observed |
|------|------------------|-------------------|
| L² | O(h³) | ~3.0 |
| H¹ | O(h²) | ~2.0 |

### P2-P1 Taylor-Hood (Stokes)

| Norm | Theoretical Order |
|------|------------------|
| L² velocity | O(h³) |
| H¹ velocity | O(h²) |
| L² pressure | O(h²) |

### P1-P1 (Stokes, unstable)

The P1-P1 pair does **not** satisfy the inf-sup condition, leading to spurious pressure oscillations and potentially incorrect velocity fields. This is included for pedagogical purposes.

---

## Key Concepts Demonstrated

1. **Saddle-point systems:** The Stokes equations lead to an indefinite (saddle-point) linear system, fundamentally different from the positive-definite systems arising from Poisson problems.
2. **Inf-sup stability:** The choice of finite element pair (P2-P1 vs P1-P1) determines whether the discrete problem is well-posed.
3. **Pseudo-elimination:** A technique for imposing Dirichlet boundary conditions that preserves matrix symmetry and sparsity.
4. **Convergence analysis:** Systematic study of numerical error as a function of mesh size to verify theoretical predictions.
5. **Flow rate computation:** Numerical integration (trapezoidal rule) to compute physically meaningful quantities from the FEM solution.
