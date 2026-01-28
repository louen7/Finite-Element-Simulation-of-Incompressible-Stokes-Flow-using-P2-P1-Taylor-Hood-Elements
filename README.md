Numerical Simulation of Incompressible Fluid Flow (Stokes P2-P1)
Project Overview
This repository contains a Finite Element Method (FEM) solver developed to simulate steady-state, viscous, incompressible fluid flows. The project focuses on solving the Stokes equations, which model fluid behavior at low Reynolds numbers by neglecting convective terms in favor of diffusion:

-nu*Delta(u) + grad(p) = f in Omega 
div(u) = 0 in Omega 

Numerical Methodology
1. Taylor-Hood (P2-P1) Elements
To ensure numerical stability and satisfy the inf-sup (Brezzi-Babuška) condition, we implement the Taylor-Hood element pair:

Velocity (u): Discretized using P2 quadratic elements associated with triangle vertices and edge midpoints.
Pressure (p): Discretized using P1 linear elements associated only with triangle vertices.

2. Block Matrix Assembly

The coupled system is assembled into a global saddle-point matrix structured by blocks:
Block [nuK, 0; 0, nuK]: Represents the diffusion of velocity components (Stiffness).
Blocks [E, F] and [G, H]: Represent the velocity-pressure coupling (Gradient and Divergence).
Quadrature: High-order Gauss-Legendre integration is used (6 points for mass, 3 points for stiffness).

3. Boundary Conditions
The solver handles both Neumann and Dirichlet conditions. We utilize the pseudo-elimination technique to impose prescribed velocity profiles at the boundaries without altering the global matrix structure.

Validation
The code's accuracy is validated through:

Poiseuille Flow: Comparison against the analytical solution for flow in a tube.

Convergence Analysis: Error estimation in L2 and H1 norms to verify theoretical orders of convergence.

Visualization: Velocity vector fields are plotted using the quiver function.

Repository Structure
principal_stokes.m: Main script for the Stokes flow simulation.
principal_dirichlet_p2 : Dirichlet Problem
principal_neumann_p2 : 
matK_elem_p2.m / matM_elem_p2.m: Computation of quadratic elemental matrices.
matE_elem.m / matF_elem.m: Computation of velocity-pressure coupling blocks.
elimine_stokes.m: Specialized routine for boundary condition enforcement.
