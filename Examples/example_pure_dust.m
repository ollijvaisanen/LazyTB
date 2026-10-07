%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% An example script for the NLTB solver involving a setup of dark matter
% only.
%
%     [1]: LazyTB: Extended LTB-solutions with Multiple Interacting Fluids.
%                  Kimmo Kainulainen, Enrico Schiappacasse, Linda Tenhu, Olli Väisänen,
%                  XX.10.2026, arXiv: XXXX.XXXXX.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Some parameters are necessary for defining the initial profiles etc.

h = 0.71;              % Hubble rate.


rb = 5*h; % Radius of the perturbation in units of Mpc/h

radius_scaling = 3.336e-4; % Unit conversion Mpc/h -> c/H0.
trb = rb*radius_scaling;

scale = 8;          % Set the grid boundary scale times the radius of the perturbation.
trgd  = trb*scale;  % Grid boundary in units of c/H0;


% -----------------------------------------------------------------------
% Define the initial density profile ------------------------------------
% -----------------------------------------------------------------------

Deltatr = 0.3*trb;
gr = @(tr) (1 - tanh((tr-trb)/Deltatr))/(1 + tanh(trb/Deltatr));

dgr = @(tr) - (cosh((tr-trb)/Deltatr)).^(-2)./(1 + tanh(trb/Deltatr))/Deltatr;
fr  = @(tr) gr(tr) + (tr/3).*dgr(tr);


% ---------------------------
% Cold dark matter profile --
% ---------------------------

% Note that while the CDM is always present in the code as a "reference
% fluid", it is not computed dynamically and does not take resources to
% run. Nothing prevents the user from setting the CDM-density fraction to
% zero, as long as the initial conditions integrals remain convergent.

Omega_CDM = 1.0;

delta_init_CDM = -1.5e-3;

gr_CDM = @(tr) delta_init_CDM * gr(tr);
fr_CDM = @(tr) delta_init_CDM * fr(tr);


% ------------------------------------------------------------------------
% Create solver instance -------------------------------------------------
% ------------------------------------------------------------------------

zin = 1000;
zout = 0;

solver = NLTBSolver(h=h,...
                    trgd = trgd,...
                    zin=zin,...
                    zout=zout,...
                    nr=350,...
                    ny=100, ...
                    Omega_CDM=Omega_CDM,...
                    gr_CDM=gr_CDM,...
                    fr_CDM=fr_CDM,...
                    Rtol=1e-8,...
                    Atol=1e-8 ...
                    );


% -----------------------------------------------------------------------
% Run the solver --------------------------------------------------------
% -----------------------------------------------------------------------

output = solver.run();

