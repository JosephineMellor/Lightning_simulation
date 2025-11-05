# Lightning Simulation in C++
Using magnetohydrodynamics with a cylindrical coordinate system to model lightning.

## Cylindrical Compressible Euler equations
Converted 1d Slope Limiting scheme (SLIC) into axi-symmetric cylindrical coordinates using source terms.
<p align="center">
  <img src="C++_scripts/euler_density.gif" alt="Sod Density Animation" width="450" />
</p>
<p align="center"><em>Figure 1: Density over Time under Cylindrical Sod Shock Test</em></p>

## Working Lighting Model (Tabulated Equation of State)
Finite difference and finite volume schemes modelling lightning strike.
- Electromagnetic source terms with subcycling.
- Lighting modelled as a cylindrically symmetric rod.
- Bilinear interpolation on a tabulated equation of state.
- Reflective Boundaries (representing the line of symmetry at r=0) at LHS of mesh and transmissive boundaries (representing the artificial open boundary with air) on RHS of mesh.
  <p align="center">
  <img src="C++_scripts/density.png" alt="Lightning Density" width="300" />
  <img src="C++_scripts/pressure.png" alt="Lightning Pressure" width="300" />
  <img src="C++_scripts/temperature.png" alt="Lightning Temperature" width="300" />
</p>
