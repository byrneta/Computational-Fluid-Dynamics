# Computational Fluid Dynamics

Here is a collection of MATLAB code (compatible with Octave) that might be of some help in solving various types of computational fluid dynamics 
problems. Functionally the codes produce valid results; however, I am sure there is room for improvement from an efficiency standpoint. The original project description from my professor is also posted for each type of problem.

### Octave Compatibility

Tested with Octave 11.3.0 (arm64 Mac), M1 Max 32GB

## Diffusion PDE

Finite difference approximation of a given couette flow between two parallel plates. 
The fluid has a constant kinematic viscosity and density. The upper plate is stationary and the lower one is suddenly set in motion with a constant velocity. Governing partial differential equation (PDE) is discretized using a first-order forward-time and second-order central space (FTCS) scheme.

See [Description](https://raw.github.com/byrneta/Computational-Fluid-Dynamics/master/diffusion/description.pdf)

### Comparative Times

Analytical = 0.05328s\
Approximation = 0.07509s

### Example Plot
![Diffusion PDE Plot](https://raw.githubusercontent.com/byrneta/Computational-Fluid-Dynamics/master/diffusion/diffusion.png)

## Convection-Diffusion PDE

Comparison between finite difference and finite volume approximations of wave propagation inside a one-dimensional channel. Fluid velocity and diffusion coefficient are given in addition to initial conditions along the channel and boundary conditions at the inlet and outlet. Forward Time-Centered Space (FTCS) and first-order upwind are used for finite volume approximations while FTCS, first-order upwind, Lax-Wendroff and MacCormack are used for finite difference.

See [Description](https://raw.github.com/byrneta/Computational-Fluid-Dynamics/master/convection-diffusion/description.pdf)

### Comparative Times

Analytical = 0.000005s

#### Finite Difference

FTCS/FTCS = 1.0439s\
Upwind/FTCS = 1.0252s\
Lax-Wendroff/FTCS = 1.5129s\
MacCormack/FTCS = 2.1235s

#### Finite Volume

FTCS/FTCS = 1.0541s\
Upwind/FTCS = 1.0274s

### Example Plots
![Convection-Diffusion PDE Plot #1](https://raw.githubusercontent.com/byrneta/Computational-Fluid-Dynamics/master/convection-diffusion/convection-diffusion-fig1.png)
![Convection-Diffusion PDE Plot #2](https://raw.githubusercontent.com/byrneta/Computational-Fluid-Dynamics/master/convection-diffusion/convection-diffusion-fig2.png)

## Elliptic PDE

Steady-state temperature distribution of a two-dimensional rectangular plate is approximated using finite difference method. Plate dimensions and boundary conditions at the edges are given. Different types of relaxation are applied: Point Successive Over-Relaxation (PSOR), Line Successive Over-Relaxation (LSOR), and Alternative Direction Implicit (ADI).

### Comparative Times

Analytical = 0.3730s\
PSOR = 4.4366s\
LSOR = 0.4166s\
ADI = 4.6181s

See [Description](https://raw.github.com/byrneta/Computational-Fluid-Dynamics/master/elliptic/description.pdf)

### Example Plots
![Elliptic PDE Plot #1](https://raw.githubusercontent.com/byrneta/Computational-Fluid-Dynamics/master/elliptic/elliptic-fig1.png)
![Elliptic PDE Plot #2](https://raw.githubusercontent.com/byrneta/Computational-Fluid-Dynamics/master/elliptic/elliptic-fig2.png)
![Elliptic PDE Plot #3](https://raw.githubusercontent.com/byrneta/Computational-Fluid-Dynamics/master/elliptic/elliptic-fig3.png)
![Elliptic PDE Plot #4](https://raw.githubusercontent.com/byrneta/Computational-Fluid-Dynamics/master/elliptic/elliptic-fig4.png)

## Vorticity-Stream Function Method

The steady-state u-velocity profile of an incompressible laminar flow within a plane channel is approximated using finite difference. The flow is governed by both vorticity and stream-function transport equations.

See [Description](https://raw.github.com/byrneta/Computational-Fluid-Dynamics/master/vorticity-streamfunction/description.pdf)

### Comparative Times

Analytical = 0.0027s\
Approximation = 292.11s

### Example Plots
![Vorticity PDE Plot #1](https://raw.githubusercontent.com/byrneta/Computational-Fluid-Dynamics/master/vorticity-streamfunction/vortstream-fig1.png)
![Vorticity PDE Plot #2](https://raw.githubusercontent.com/byrneta/Computational-Fluid-Dynamics/master/vorticity-streamfunction/vortstream-fig2.png)
![Vorticity PDE Plot #3](https://raw.githubusercontent.com/byrneta/Computational-Fluid-Dynamics/master/vorticity-streamfunction/vortstream-fig3.png)
![Vorticity PDE Plot #4](https://raw.githubusercontent.com/byrneta/Computational-Fluid-Dynamics/master/vorticity-streamfunction/vortstream-fig4.png)
![Vorticity PDE Plot #5](https://raw.githubusercontent.com/byrneta/Computational-Fluid-Dynamics/master/vorticity-streamfunction/vortstream-fig5.png)
![Vorticity PDE Plot #6](https://raw.githubusercontent.com/byrneta/Computational-Fluid-Dynamics/master/vorticity-streamfunction/vortstream-fig6.png)
