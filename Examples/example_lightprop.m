%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% An example script for the light propagation solver
%
%     [1]: LazyTB: Extended LTB-solutions with Multiple Interacting Fluids.
%                  Kimmo Kainulainen, Enrico Schiappacasse, Linda Tenhu, Olli Väisänen,
%                  XX.10.2026, arXiv: XXXX.XXXXX.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Either run the NLTB-solver or load a pre-calculated metric solution from
% a file. "output" is a struct with three fields, "solution", "derived" and "input". 
% See the NLTB-example for more details.

load("output.mat","output"); % Just an example, add your own here.

% ------------------------------------------------------------------------
% Metric parameters contained in a (Ny,Nr)-arrays. ------------------------
% ------------------------------------------------------------------------

% Grid parameters
metric.y  = output.solution.y;    % y = ln(a_bg), used as a time variable. 
                                  % This should extend by at least 2*da beyond the observer location
metric.tr = output.solution.tr;   % tr = radius in units of c/H0

% Metric components. See appendix A.1 in [1] for definitions.
metric.detat = output.solution.detat; 
metric.e     = output.solution.e;
metric.dbt   = output.solution.dbt;
metric.b     = output.derived.b;


% The background fluid parameters. This choice corresponds to a model with
% CDM, radiation and a cosmological constant.
omega_M = 0.3;
omega_Rad = 4.8993e-05;
omega_Lambda = 1 - omega_M - omega_Rad;

                 
% Background adot and adotdot in units of H0 and H0^2.
metric.adot    = @(a) sqrt(omega_M*a.^(-1) + omega_Rad*a.^(-2) + omega_Lambda*a.^2 ); 
metric.adotdot = @(a) -( omega_M*a.^(-2) + 2*omega_Rad*a.^(-3) - 2*omega_Lambda*a)/2;

% ------------------------------------------------------------------------
% Observer parameters ----------------------------------------------------
% ------------------------------------------------------------------------

% Observer background scale factor. NOTE: Due to central differences this
% should be smaller than the final scale factor of the NLTB-solution minus
% 2*da. If you want to run the light propagation for ainit = 1.0, you need
% to run the NLTB-solution a bit to the future.
ainit = 1.0;
aend = 1/(1 + 10);


% Observer distance from the origin in units of 1/H0.
rb = metric.tr(end)/8;
r0 = 0.5 * rb; 
%r0 = 5.0 * rb; 
%r0 = 3.0 * rb;


% Steps for numerical derivatives inside the solver. Note that da is NOT
% directly the timestep for the ode-solver. The solver uses ode45 with
% adaptive timestep instead.
da = 1e-3;
dr = output.solution.tr(2) - output.solution.tr(1); 

% Solver is run for Nrays equally spaced rays between minimum and maximum
% angles. Angle = 0 is directly towards the origin.
% NOTE: The calculation is performed in spherical coordinates and setting
% angle = 0 will run into the coordinate singularity. This can be avoided
% by including a small non-zero angle. Passing close to the origin does not
% appear to add significant numerical error, as long as the metric is
% regular.
angles = linspace(1e-5,pi,200);


% ---------------------------------------------------------------------
% Create the solver instance -------------------------------------------
% ----------------------------------------------------------------------

LP_solver = LightPropSolver(metric,...
                            ainit = ainit,...
                            aend = aend,...
                            r0 = r0,...
                            da = da,...
                            dr = dr,...
                            angles=angles);

% ------------------------------------------------------------------------
% Run the solver ---------------------------------------------------------
% ------------------------------------------------------------------------

LP_solution = LP_solver.run();

% Compute differences in solution to the background FLRW-model.
LP_diff = LP_solver.get_output_differences(LP_solution);
