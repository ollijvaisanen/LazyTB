%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% An example script for the NLTB solver involving a setup of interacting 
% dark matter and dark energy.
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


% ---------------------------
% Cold dark matter profile --
% ---------------------------

% Note that while the CDM is always present in the code as a "reference
% fluid", it is not computed dynamically and does not take resources to
% run. Nothing prevents the user from setting the CDM-density fraction to
% zero, as long as the initial conditions integrals remain convergent.

% Set the non-interacting CDM fraction to zero.
Omega_CDM = 0.0;

delta_init_CDM = -1.5e-3;

gr_CDM = @(tr) delta_init_CDM * gr(tr);
fr_CDM = @(tr) delta_init_CDM * fr(tr);


% -----------------------------------------------------------------------
% New interacting dark matter fluid -------------------------------------
% -----------------------------------------------------------------------

% We cannot add interactions to the "reference CDM-fluid". Instead, we create
% a new fluid with w=0.

Omega_DM = 0.3;

gr_DM = @(tr) delta_init_CDM * gr(tr);
fr_DM = @(tr) delta_init_CDM * fr(tr);


% -----------------------------------------------------------------------
% Dynamical DE ----------------------------------------------------------
% -----------------------------------------------------------------------

Omega_DE = 0.7;

alpha = -1.15e-5;  
wba = -1 + 1e-5;

% This time, we set the initial DE perturbation to zero.
delta_init_DE = 0;

gr_DE = @(r) delta_init_DE/delta_init_CDM * gr_CDM(r);
fr_DE = @(r) delta_init_DE/delta_init_CDM * fr_CDM(r);

% -----------------------------------------------------------------------
% Set the non-CDM fluid parameters ----------------------------------------
% -----------------------------------------------------------------------

Omegas   = [Omega_DM,Omega_DE];      % The z=0 -energy density fractions of each non-CDM fluid
fluid_names = {"dark_matter","dark_energy"};  % The names of each non-CDM fluid. These will only be used in the output struct.

% The initial profiles of each non-CDM fluid, in the same order as above.
grs = {gr_DM,gr_DE};  
frs = {fr_DM,fr_DE};

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
    
    ny = sz(1);
    nr = sz(2);

    y = reshape(log(a),[ny,1]);

    cs2_DM  = zeros([1,sz]);
    w_DM    = zeros([1,sz]);
    wdot_DM = zeros([1,sz]);
    wr_DM   = zeros([1,sz]);
    wbg_DM  = zeros([1,sz]);
    wint_DM = reshape(repmat(y,[1,nr]),[1,sz]);

    Delta_DE = reshape(Deltas(3,:),[1,sz]); % DE is the 3rd fluid, after the reference CDM and interacting DM.
    w_DE     = wba*Delta_DE.^(alpha-1);
    cs2_DE   = alpha*w_DE;
    wdot_DE  = 3*(1+wba)*(alpha-1)*w_DE;
    wr_DE    = 0*w_DE;
    wbg_DE   = wba.*ones([1,sz]);
    wint_DE  = reshape(repmat((1 + wba)*y,[1,nr]),[1,sz]);


    cs2  = [cs2_DM ; cs2_DE];
    w    = [w_DM   ; w_DE];
    wdot = [wdot_DM; wdot_DE];
    wr   = [wr_DM  ; wr_DE];
    wbg  = [wbg_DM ; wbg_DE];
    wint = [wint_DM; wint_DE];
end
EOSPMS = @(a,Deltas,sz) EOSPMS_DE(a,Deltas,sz,wba,alpha);


% ------------------------------------------------------------------------
% Define interaction terms -----------------------------------------------
% ------------------------------------------------------------------------


function [Qt,arQr] = interactions(a,v,Delta,w,cs2)
    zeta_scaled = 1;
    signs = [1;-1];


    Qt = zeta_scaled * signs * ( 1/sqrt(1 - v(3)^2) - 1/sqrt(1 - v(2)^2) );
    arQr = zeta_scaled * signs * ( 1/sqrt(1 - v(3)^2) * v(3) - 1/sqrt(1 - v(2)^2) * v(2) );
    
end

% ----------------------------------------------------------------------
% Numerical diffusion --------------------------------------------------
% ----------------------------------------------------------------------

% Numerical diffusion becomes necessary when the sound horizon of some
% fluid becomes large compared with the size of the grid. Note that even
% when the diffusion has no visible effect on the output solution, a small
% amount of diffusion can have a drastic effect on the computation time.

                          
DDAll = 10^-7;                            % Diff reg scale for the fluid-2 (log) density
DvAll = 10^-7;                            % Diff reg scale for the fluid-2 velocity 
 

Dorigin = 1; Dsmall = 1.0; width = 0.02;

famp  = @(z) (1/3)/sqrt(0.3*(1+z)^3);

% Diffusion has to be set to zero for the DM-fluid, as the DM distribution 
% is at rest apart from the small interactions. As such, any diffusion 
% would result in a dominant term. The effects of adding diffusion here
% are easy to see by comparing the resulting DM density contrast to the
% analytical CDM density contrast computed from the metric (in output.derived).
DDFlav = DDAll.*[0;1.0];
DvFlav = DvAll.*[0;1.0];

DDiff = @(z,tr)  DDFlav*famp(z)*( (tr<2*trb).*( Dorigin*exp(-(tr/(width*2*trb)).^2) + Dsmall  ) + Dsmall*(tr>=2*trb).*(tr/(2*trb)).^5 );
vDiff = @(z,tr)  DvFlav*famp(z)*( (tr<2*trb)*(tr/(2*trb)).^5 + (tr>=2*trb).*(tr/(2*trb)).^5 );

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
                    Omegas = Omegas,...
                    fluid_names = fluid_names,...
                    EOSPMS = EOSPMS,...
                    grs = grs,...
                    frs = frs,...
                    Rtol=1e-8,...
                    Atol=1e-8, ...
                    DDiff = DDiff,...
                    vDiff = vDiff,...
                    interactions = @interactions...
                    );


% -----------------------------------------------------------------------
% Run the solver --------------------------------------------------------
% -----------------------------------------------------------------------

output = solver.run();
