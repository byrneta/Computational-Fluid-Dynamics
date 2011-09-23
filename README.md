# Computational Fluid Dynamics

Here is a collection of MATLAB code that might be of some help in solving various types of computational fluid dynamics 
problems. I wrote them during graduate school while pursuing an M.S. in Mechanical Engineering and functionally the codes produce valid results; however, I am sure there is room for improvement from an efficiency standpoint. The original project description from my professor is also posted for each type of problem.

## Diffusion PDE

Finite difference approximation of a given couette flow between two parallel plates. 
The fluid has a constant kinematic viscosity and density. The upper plate is stationary and the lower one is suddenly set in motion with a constant velocity. Governing PDE is discretized using a first-order forward-time and second-order central space (FTCS) scheme.

See [Description](https://raw.github.com/byrneta/Computational-Fluid-Dynamics/master/diffusion/description.pdf)

## Convection-Diffusion PDE

Comparison between finite difference and finite volume approximations of wave propagation inside a one-dimensional channel. Fluid velocity and diffusion coefficient are given in addition to initial conditions along the channel and boundary conditions at the inlet and outlet. FTCS and first-order upwind are used for finite volume approximations while FTCS, first-order upwind, Lax-Wendroff and MacCormack are used for finite difference.

See [Description](https://raw.github.com/byrneta/Computational-Fluid-Dynamics/master/convection-diffusion/description.pdf)

## Elliptic PDE

Steady-state temperature distribution of a two-dimensional rectangular plate is approximated using finite difference method. Plate dimensions and boundary conditions at the edges are given. Different types of relaxation are applied: PSOR, LSOR, and ADI.

See [Description](https://raw.github.com/byrneta/Computational-Fluid-Dynamics/master/elliptic/description.pdf)

## Vorticity-Stream Function Method

The steady-state u-velocity profile of an incompressible laminar flow within a plane channel is approximated using finite difference. The flow is governed by both vorticity and stream-function transport equations.

See [Description](https://raw.github.com/byrneta/Computational-Fluid-Dynamics/master/vorticity-streamfunction/description.pdf)
