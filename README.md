# LazyTB_testing

IN PROGRESS: If you somehow end up here, the code will appear tomorrow.

LazyTB is publicly available MATLAB code for computing cosmological observables in spherically symmetric spacetimes with an arbitrary number of ideal fluids. The code comprisses of a solver for the Einstein field equations and the fluid equations, routines for specifying the initial conditions and diagnostics tools for the numerical error. In addition, the code includes a solver for observables such as angular diameter distances along light rays for an observer at an arbitrary location of the spacetime.

A python version of the code is under development and a link to the repository will be added here.

More information on the details of the code can be found in. 

[1] LazyTB: Extended LTB-solutions with Multiple Interacting Fluids. Kimmo Kainulainen, Enrico Schiappacasse, Linda Tenhu, Olli Väisänen, 7.10.2026, arXiv: XXXX.XXXXX.

Please cite this article when using our code.

# Quick start
Install Matlab 2021a or newer, (the code was built and tested on Matlab2024a). Download the contents in /LazyTB to your local folder or add their location to your MATLAB path.

# NLTBSolver

Note that the calculation is performed and the output is given in terms of perturbations from a background FLRW solution, whose energy components are defined by the density fractions Omega_CDM and Omegas, and whose equations of state (EOS) are given in the EOSPMS-function, (see the definitions below).

### Arguments to NLTBSolver:

##### h
Scalar. Background Hubble rate is $H_0 = h\cdot$ 100 Mpc/km/s.

##### trgd
Scalar. The size of the radial grid in units of 1/H0.

##### Omega_CDM
Scalar. Density fraction of a pressureless fluid. The computation is performed in the rest frame of this fluid and its evolution is handled analytically. Omega_CDM can be set to zero, as long as the initial conditions integrals (see ch. 2.5 of [1]) are convergent.

##### gr_CDM
Function handle $$g(r)$$, $r$ is in units of $H_0^{-1}$. Radial density profile of the pressureless fluid. See equation (2.38) and (2.42) of [1].

##### fr_CDM
Function handle $$f(r)$$, $r$ is in units of $H_0^{-1}$. Radial density profile of the pressureless fluid. See equation (2.40) of [1]. This function has to equal $$f(r) = g(r) + (r/3)g'(r)$$.

##### zin 
Scalar. Initial redshift

##### zout
Scalar. Final redshift

##### nr, ny
Scalar. The number of grid points in the radial and time direction in the output array

##### Rtol, Atol 
Scalar. Relative and absolute error tolerance passed to the pdepe-solver.

##### DDiff, vDiff
Function handles $D_i(z,r)$, where $r$ is in units of $H_0^{-1}$. Diffusion functions for density and velocity perturbations. See equations (??). These should be functions DDiff(z,r) and the output shape must
be either scalar or (Nfluids,1).

### Adding fluid components
LazyTB supports adding an arbitrary number of additional fluid components with user-defined equations of state and interaction terms. An analytically-handled pressureless CDM-fluid is always present, although its density can be set to zero as long as the initial conditions converge (see Ch. 2.5 in [1]). In the below definitions, Nfluids is the number of non-CDM fluid components.

##### fluid_names
Cell array of strings with size (Nfluids,1). The name of each non-CDM fluid in the output struct

##### Omegas
Array with size (Nfluids,1). The background density fraction of each non-CDM fluid.

##### grs, frs
Cell arrays containing function handles. See gr_CDM and fr_CDM above.

##### EOSPMS
A function handle returning the equations-of-state parameters of the fluid. The arguments must be as follows, where the input size (Ny,Nr) is either (1,1) or (ny,nr).

| Argument | Definition                                                                                                |
|----------|-----------------------------------------------------------------------------------------------------------|
|    a     | Background scale factor, used as a time variable                                                          |
| Deltas   | (Nfluids+1,Ny,Nr)-array containing the density contrast of each fluid at a given time and radius. |
|   sz     | Array [Ny,Nr]. This argument only exists to get around matlab squeezing trailing dimensions.              |

The outputs must all be of the size (Nfluids,Ny,Nr). Note that if there are no non-CDM fluids, the outputs should be empty arrays and if there is only a single non-CDM fluid, the leading singleton dimension should be kept.

| Output   | Definition                                                                                             |
|----------|--------------------------------------------------------------------------------------------------------|
|   cs2    | The speed of sound squared, defined as in eq. (2.29) of [1].                                           |
|   w      | The EOS parameter $p = w\rho$.                                                                         |
|   wdot   | The partial derivative of $w$ wrt. coordinate time                                                     |
|   wr     | The partial radial derivative.                                                                         |
|   wbg    | The background value of $w$. Our boundary conditions require $w_{bg}\to w$ at the outer grid boundary. |
|   wint   | Value of integral defined below at scale factor a.                                                 |

The last output of EOSPMS "wint" should return the values of the integral
```math
w_{int} = \int_0^y (1 + w_{bg}(y'))dy' 
```
where $y = ln(a)$. This is used to compute the evolution of background energy densities. See eq. (A.4) of [1].

##### Interactions

Function handle describing the interactions. The inputs should be as follows:

| Field    | Definition                                                                                            |
|----------|-------------------------------------------------------------------------------------------------------|
|    a     | Background scale factor, used as a time variable                                                      |
| Deltas   | (Nfluids+1,1)-array containing the density contrast of each non-CDM fluid at a given time and radius. |
|    w     | (Nfluids+1,1)-array containing the EOS parameter $p = w \rho$ for each fluid.                         |
|   cs2    | (Nfluids+1,1)-array containing the speed of sound of each fluid                                       |

The output must be two (Nfluids,1) arrays containing the interaction vectors $\bar Q^t$ and $a_r \bar Q^r$, see eqs (A.8-A.9) in [1].


### Running of the NLTBSolver
The NLTBSolver is run with the method solver.run(add_input=true,add_derived=true). The keyword argument specifies whether the solver saves additional information necessary for the light propagation solver. The output has the following structure:

#### output.solution
Array containing the output. See appendix (A.1) of [1] for all definitions. The output fields are as follows:

| Field                | Definition                                                                                               |
|----------------------|----------------------------------------------------------------------------------------------------------|
| detat                | Perturbation in dimensionless angular metric component $\eta_\theta = 1 + \delta\eta_\theta$.            |
| e                    | Dimensionless curvature $e$.                                                                             |
| dbt                  | Perturbation in local expansion rate $b_\theta = -K^\theta_{\theta}/\bar H = 1 + \delta b_\theta$        |
| logDelta_fluid_name  | Logarithm of density contrast. One field for each fluid                                                  |
| v_fluid_name         | Peculiar velocity of each fluid.                                                                         |
| tr                   | Radial grid in units of $H_0^{-1}$.                                                                      |
| y                    | Time grid expressed in ln(a).                                                                            |

#### output.input
A copy of the input struct. Only if add_input=true.

#### output.derived
A struct for derived quantities. See appendix (A.1-A.2) of [1] for all definitions. Only if add_derived=true.

| Field                | Definition                                                                                               |
|----------------------|----------------------------------------------------------------------------------------------------------|
|      br_mom          | The dimensionless extrinsic curvature component $b_r$ computed using the momentum constraint.            |
|      br_ham          | The dimensionless extrinsic curvature component $b_r$ computed using the hamiltonian constraint.         |
|     Delta_CDM        | The density contrast of pressureless CDM, computed analytically using eq. (A.13) of [1].                 |
|       hb             | The background FLRW expansion rate in units of $H_0$.                                                    |
|        b             | The local expansion rate $b = b_r + 2b_\theta$, computed using the momentum constraint.                  |


# LightPropSolver

### Arguments of LightPropSolver

##### metric
A struct containing the metric for the computation in spherical coordinates. The input is of the same form as the output of NLTBSolver, but it is not necessary. The metric struct must have the following fields:

| Field                | Definition                                                                                                       |
|----------------------|------------------------------------------------------------------------------------------------------------------|
| y                    | Time grid for metric variables expressed in ln(a).                                                               |
| tr                   | Radial grid for metric variables in units of $H_0^{-1}$.                                                         |
| detat                | Perturbation in dimensionless angular metric component $\eta_\theta = 1 + \delta\eta_\theta$.                    |
| e                    | Dimensionless curvature $e$.                                                                                     |
| dbt                  | Perturbation in local expansion rate $b_\theta = -K^\theta_{\theta}/\bar H = 1 + \delta b_\theta$                |
| b                    | The local expansion rate $b = b_r + 2b_\theta$, computed using the momentum constraint.                          |
| adot                 | A function handle $\dot a(y)$ returning the derivative of the background scale factor in units of $H_0$.         |
| adotdot              | A function handle $\ddot a(y)$ returning the 2nd derivative of the background scale factor in units of $H_0^2$.  |

##### ainit
Scalar. The scale factor of the observer. Note that the calculation is done backward in time and due to the use of central derivatives inside the solver, the metric must extend by 2*da beyond the observer timeslice.

##### aend
Scalar. The end scale factor

##### r0, phi0
Scalar. The observer radius and angle.


##### da, dr
The steps used in computing finite difference derivatives of the metric inside the solver. Note that the optical equations themselves are solved with an adaptive solver and the timestep used here is not directly related to it.

##### angles
Angles of the outgoing light rays from the observer. This can be an array for solving a set of rays simultaneously. Zero angle is directed towards the origin.

##### ode_settings
ODE settings passed on to ode45 in the standard form.

### Running the LightPropSolver
The solver is run with the method LightPropSolver.Run(). The output struct contains the following fields:

| Field                | Definition                                                                                               |
|----------------------|----------------------------------------------------------------------------------------------------------|
| a                    | The background scale factor at the solution timesteps.                                                   |
| rrays                | The radius of each ray at the different timesteps in units of $H_0^{-1}$.                                |
| phirays              | The angle of each ray at the different timesteps. Zero angle points towards the origin.                  |
| s2                   | The components of the $s_2$ frame vector. See Ch. 2.6 in [1].                                            |
| expansion            | The expansion scalar. See Ch. 2.6 in [1].                                                                |
| angular_distance     | The angular diameter distance.                                                                           |
| k_phi                | The angular component of the momentum vector.                                                            |
| zrays                | The redshift observed along each ray.                                                                    |
| k0                   | The timelike component of the momentum vector.                                                           |
| angles               | The angles-vector appended on the solution for convenience.                                              |

Differences from the quantities observed in the background FLRW-spacetime can be obtained by running the method LightPropSolver.get_output_differences(LP_solution). The output has the following fields:

| Field                     | Definition                                                                                          |
|---------------------------|-----------------------------------------------------------------------------------------------------|
| zrays_relative            | Relative difference in redshift along each ray.                                                     |
| temperature_relative      | Relative difference in observed radiation temperature along each ray.                               |
| angular_distance_relative | Relative difference in angular distance along each ray.                                             |
| angular_distance_bg       | Angular distance observed in the background FLRW-spacetime.                                         |


# Examples

The scripts necessary for running the cases of Ch. 3 of [1] are contained in the Examples-folder. These examples highlight the following parts of the code:

1) example_pure_dust: Running the NLTBSolver without any additional fluids.
2) example_dynamical_DE: Adding a non-trivial fluid to the solver.
3) example_DM_rad: Adding a large-sound speed fluid to the solver with numerical diffusion. This example should be run with a larger precision and it will take quite a bit longer to run than the previous two cases.
4) example_DM_rad_lambda: Adding multiple fluids to the EOSPMS function. As before, this run takes a bit longer due to the large sound speeds.
5) example_interacting_DM_DE: Adding interactions between the fluids.
6) example_lightprop: Running the LightPropSolver, save one of the solutions in the previous examples for use here.


