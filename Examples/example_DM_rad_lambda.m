%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% An example script for the NLTB solver involving a setup of dark matter,
% radiation and a cosmological constant.
%
%     [1]: LazyTB: Extended LTB-solutions with Multiple Interacting Fluids.
%                  Kimmo Kainulainen, Enrico Schiappacasse, Linda Tenhu, Olli Väisänen,
%                  XX.10.2026, arXiv: XXXX.XXXXX.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Some parameters are necessary for defining the initial profiles etc.
h = 0.71;              % Hubble rate.

rb = 50*h; % Radius of the perturbation in units of Mpc/h

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


% -----------------------------------------------------------------------
% Cold dark matter profile ----------------------------------------------
% -----------------------------------------------------------------------

% Note that while the CDM is always present in the code as a "reference
% fluid", it is not computed dynamically and does not take resources to
% run. Nothing prevents the user from setting the CDM-density fraction to
% zero, as long as the initial conditions integrals remain convergent.

Omega_CDM = 0.3;

delta_init_CDM = -1e-6;

gr_CDM = @(tr) delta_init_CDM * gr(tr);
fr_CDM = @(tr) delta_init_CDM * fr(tr);



% ------------------------------------------------------------------------
% Radiation --------------------------------------------------------------
% ------------------------------------------------------------------------

% Set the radiation energy density.
Trad = 2.725;
rhoc   = 8.0992e-47*h^2;           % Critical density
rhog   = 2.0747e-51*(Trad/2.75)^4; % Radiation energy density as a function of temperature
Omega_rad = rhog/rhoc;

delta_init_rad = 4/3 * delta_init_CDM;
gr_rad = @(tr) delta_init_rad * gr(tr);
fr_rad = @(tr) delta_init_rad * fr(tr);


% -----------------------------------------------------------------------
% Cosmological constant -------------------------------------------------
% -----------------------------------------------------------------------

Omega_Lambda = 1 - Omega_CDM - Omega_rad;


gr_Lambda = @(tr) zeros(size(tr));
fr_Lambda = @(tr) zeros(size(tr));


% -----------------------------------------------------------------------
% Set the other fluid parameters ----------------------------------------
% -----------------------------------------------------------------------

Omegas   = [Omega_rad,Omega_Lambda];      % The z=0 -energy density fractions of each non-CDM fluid
fluid_names = {"photons","dark_energy"};  % The names of each non-CDM fluid. These will only be used in the output struct.

% The initial profiles of each non-CDM fluid, in the same order as above.
grs = {gr_rad,gr_Lambda};  
frs = {fr_rad,fr_Lambda};

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
function [cs2,w,wdot,wr,wbg,wint] = EOSPMS_DM_rad_Lambda(a,Deltas,sz)

    cs2  = repmat([1/3;0] ,[1,sz]);
    w    = repmat([1/3;-1],[1,sz]);
    wdot = zeros([2,sz]);
    wr   = zeros([2,sz]);
    wbg  = repmat([1/3;-1],[1,sz]);

    Na = sz(1);
    Nr = sz(2);
    wint_DE      = reshape(zeros(size(a)),[1,Na]);
    wint_photons = reshape(4/3 * log(a)  ,[1,Na]);
    wint         = repmat([wint_photons;wint_DE],[1,1,Nr]);
end
EOSPMS = @(a,Deltas,sz) EOSPMS_DM_rad_Lambda(a,Deltas,sz);


% ----------------------------------------------------------------------
% Numerical diffusion --------------------------------------------------
% ----------------------------------------------------------------------

% Numerical diffusion becomes necessary when the sound horizon of some
% fluid becomes large compared with the size of the grid. Note that even
% when the diffusion has no visible effect on the output solution, a small
% amount of diffusion can have a drastic effect on the computation time.

% This example needs a bit heavier diffusion to keep the radiation fluid
% in check. 
DDAll = 10^-3;                            % Diff reg scale for the fluid-2 (log) density
DvAll = 10^-3;                            % Diff reg scale for the fluid-2 velocity 

Dorigin = 0.0; Dsmall = 1.0; width = 0.02;

famp  = @(z) (1/3)/sqrt(4.9e-5*(1+z)^4 + 0.3*(1+z)^3);

DDFlav = DDAll;
DvFlav = DvAll;

DDiff = @(z,tr)  DDFlav*famp(z)*( (tr<2*trb).*( Dorigin*exp(-(tr/(width*2*trb)).^2) + Dsmall  ) + Dsmall*(tr>=2*trb).*(tr/(2*trb)).^5 );
vDiff = @(z,tr)  DvFlav*famp(z)*( (tr<2*trb)*(tr/(2*trb)).^5 + (tr>=2*trb).*(tr/(2*trb)).^5 );


% -----------------------------------------------------------------------
% Optionally freeze the radiation density perturbation at ---------------
% a > freeze_scale_factor -----------------------------------------------
% -----------------------------------------------------------------------

% Both of these parameters can be (Nfluids,1)-column vectors if the system
% contains multiple fluids.

freeze_scale_factor = inf;
freeze = 1;


% ------------------------------------------------------------------------
% Create solver instance -------------------------------------------------
% ------------------------------------------------------------------------

zin = 1e5;       % Initial redshift
zend = 1/1.1-1;  % End redshift. We set this to zend < 0 so that central derivatives function at z = 0 in the light ray calculation.


solver = NLTBSolver(h=h,...
                    trgd = trgd,...
                    zin=zin,...
                    zout=zend,...
                    nr=300,...
                    ny=150, ...
                    Omega_CDM=Omega_CDM,...
                    gr_CDM=gr_CDM,...
                    fr_CDM=fr_CDM,...
                    Omegas = Omegas,...
                    fluid_names = fluid_names,...
                    EOSPMS = EOSPMS,...
                    grs = grs,...
                    frs = frs,...
                    Rtol=1e-13,...
                    Atol=1e-13, ...
                    DDiff = DDiff,...
                    vDiff = vDiff,...
                    freeze = freeze,...
                    freeze_scale_factor = freeze_scale_factor...
                    );

% Note: We can leave any argument out if we want to use a default value. 

% -----------------------------------------------------------------------
% Run the solver --------------------------------------------------------
% -----------------------------------------------------------------------

output = solver.run();

