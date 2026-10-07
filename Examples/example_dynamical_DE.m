%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% An example script for the NLTB solver involving a setup of dark matter
% and dynamical dark energy.
%
%     [1]: LazyTB: Extended LTB-solutions with Multiple Interacting Fluids.
%                  Kimmo Kainulainen, Enrico Schiappacasse, Linda Tenhu, Olli Väisänen,
%                  XX.10.2026, arXiv: XXXX.XXXXX.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Some parameters are necessary for defining the initial profiles etc.
h = 0.71;              % Hubble rate.

rb_old = 10*h;           % Radius of the LTB-bubble in units of Mpc/h. 

radius_scaling = 3.336e-4; % Unit conversion Mpc/h -> c/H0.  
trb_old = rb_old * radius_scaling; % Convert the radius to units of c/H0

scale_old = 4; % The boundary of the grid is set at trb*scale.

%----------------------------------------------------------------------

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

Omega_CDM = 0.3;

normalization_correction = 1.0344; % Difference in normalizations between us and Marra et al.
delta_init_CDM = 1.5e-3 * normalization_correction;

gr_CDM = @(tr) delta_init_CDM * gr(tr);
fr_CDM = @(tr) delta_init_CDM * fr(tr);


% -----------------------------------------------------------------------
% Dynamical DE ----------------------------------------------------------
% -----------------------------------------------------------------------

Omega_DE = 0.7;

%alpha = 0;         % Zero sound speed
%alpha = -1.15e-5;  % cs ~ 0.003
alpha = -1.15e-3;   % cs ~ 0.03, this one needs a small diffusion term.

wba   = -0.8;

%
% The EOS parameters function. 
%
% Input:      a: Background scale factor. Scalar or (Ny,1)-array.
%        Deltas: Array of fluid density contrasts, (Nfluids,Ny,Nr)-array.
%                The first fluid is always CDM, followed by dynamical fluids
%                in the same order as in Omegas.
%            sz: Array in the size [ny,nr]. This argument only exists because
%                Matlab wants to squeeze trailing dimensions.
%
% Output: Each array should be of the shape (Nfluids,nr,ny). Note that if Nfluids = 1,
%         leave the leading singleton dimension in place!
%
%           cs2: Sound speed of each fluid. See eq. (2.29) in [1].
%             w: EOS parameter p = w*rho.
%          wdot: Time-partial derivative of w with r and rho constant.
%            wr: Radial partial derivative of w with t and rho constant
%           wbg: Background EOS parameter.
%          wint: Integral of (1 + wbg(y))dy from 0 to y(a), where y = ln(a).
%                This is used to compute background energy density of each fluid
%                component. If wbg = constant, wint = (1 + wbg)*ln(a).
%
function [cs2,w,wdot,wr,wbg,wint] = EOSPMS_DE(a,Deltas,sz,wba,alpha)
    
    Delta_DE = reshape(Deltas(2,:),[1,sz]);
    w     = wba*Delta_DE.^(alpha-1);
    cs2   = alpha*w;
    wdot  = 3*(1+wba)*(alpha-1)*w;
    wr    = 0*w;
    wbg   = wba.*ones(sz);

    wint = reshape((1 + wbg).*log(a).',[1,sz]);
end
EOSPMS = @(a,Deltas,sz) EOSPMS_DE(a,Deltas,sz,wba,alpha);


delta_init_DE = delta_init_CDM * (1 + wba) * (5 - 6*alpha*wba) / (5 - 15*wba + 9*alpha*wba);

gr_DE = @(r) delta_init_DE/delta_init_CDM * gr_CDM(r);
fr_DE = @(r) delta_init_DE/delta_init_CDM * fr_CDM(r);

% -----------------------------------------------------------------------
% Set the non-CDM fluid parameters ----------------------------------------
% -----------------------------------------------------------------------

Omegas   = [Omega_DE];          % The z=0 -energy density fractions of each non-CDM fluid
fluid_names = {"dark_energy"};  % The names of each non-CDM fluid. These will only be used in the output struct.

% The initial profiles of each non-CDM fluid, in the same order as above.
grs = {gr_DE};  
frs = {fr_DE};



% ----------------------------------------------------------------------
% Numerical diffusion --------------------------------------------------
% ----------------------------------------------------------------------

% Numerical diffusion becomes necessary when the sound horizon of some
% fluid becomes large compared with the size of the grid. Note that even
% when the diffusion has no visible effect on the output solution, a small
% amount of diffusion can have a drastic effect on the computation time.

% These values seem ok for the case "cs ~ 0.03". For the others the diffusion 
% should be set to zero.

                            
DvAll = 10^-7;                            % Cases A-C. Diff reg scale for the fluid-2 velocity 
%DvAll = 10^-6;                           % Case D

DDAll = 10^-7; Dorigin = 1e4; Dsmall = 1.0; width = 0.02;     % Case A. Diff reg scale for the fluid-2 (log) density
%DDAll = 10^-6; Dorigin = 1e3; Dsmall = 1.0; width = 0.03;    % Case B. Diff reg scale for the fluid-2 (log) density
%DDAll = 10^-5; Dorigin = 1e3; Dsmall = 1.0; width = 0.03;    % Case C-D. Diff reg scale for the fluid-2 (log) density


famp  = @(z) (1/3)/sqrt(0.3*(1+z)^3);

% If you had more than one non-CDM fluid, you could set the diffusion
% functions separately. In this example, DDFlav is just a number.
DDFlav = DDAll;
DvFlav = DvAll;

DDiff = @(z,tr)  DDFlav'*famp(z)*( (tr<2*trb).*( Dorigin*exp(-(tr/(width*2*trb)).^2) + Dsmall  ) + Dsmall*(tr>=2*trb).*(tr/(2*trb)).^5 );
vDiff = @(z,tr)  DvFlav'*famp(z)*( (tr<2*trb)*(tr/(2*trb)).^5 + (tr>=2*trb).*(tr/(2*trb)).^5 );


% ------------------------------------------------------------------------
% Create solver instance -------------------------------------------------
% ------------------------------------------------------------------------

zin  = 1000;
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
                    Omegas = Omegas,...
                    fluid_names = fluid_names,...
                    EOSPMS = EOSPMS,...
                    grs = grs,...
                    frs = frs,...
                    Rtol=1e-8,...
                    Atol=1e-8, ...
                    DDiff = DDiff,...
                    vDiff = vDiff...
                    );


% -----------------------------------------------------------------------
% Run the solver --------------------------------------------------------
% -----------------------------------------------------------------------

output = solver.run();

