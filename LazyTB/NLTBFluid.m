%---------------------------------------------------------------------------------------------------
% LTB-model wVth N fluids: 
%--------------------------------------------------------------------------------------------------- 
%
% Created: 25.10.2013
% Last modified: 17.12.2013 
%---------------------------------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Solver function for the ADM -equations of motion
% ------------------------------------------------
%
% The input structs should contain the following fields:
%
%
% input: Original input struct
%
%              h: Current Hubble rate is set to H0 = h*100 km/s/Mpc.
%           trgd: Outer grid boundary in units of c/H0.
%            zin: Initial background redshift
%           zout: End background redshift
%             nr: Number of radial grid points
%             ny: Number of temporal grid points in output
%           Rtol: Relative tolerance passed to pdepe
%           Atol: Absolute tolerance passed to pdepe
%           DDiff: Handle for the Delta-diffusion function. Input as DDiff(z,r), 
%                  output shape (Nfl,1). See Appendix A.5. in [1].
%           vDiff: Same as above for the v-diffusion function.
%               y: Grid for the time variable y = ln(a), (1,Ny)-array.
%              tr: Grid for the radial variable in units of c/H0, (1,Nr)-array
%         Nfluids: Number of non-CDM fluids
%     fluid_names: {Nfluids,1}-cell array of fluid names.
%       Omega_vec: (1,Nfluids+1)-array containing current fluid density fractions, 
%                  with CDM included.
%          gr_vec: A cell array containing initial conditions g(r)-functions each fluid. See eq. (2.38) in [1].
%          fr_vec: Same as above for the function f(r) = g(r) + (r/3)g'(r). See eq. (2.40) in [1].
%          EOSPMS: A function returning the EOS parameters for each fluid.
%    interactions: A function defining the interaction vectors between the fluids. See eq: (A.8-A.9) in [1].
% freeze_scale_factor: Optionally freeze the fluid equations at a > freeze_scale_factor. Use this to conserve computation time.
%          freeze: Set this to 0 to enable freezing and 1 otherwise. This can be a (Nfluids,1)-array.
%
%           
%
%             [1]: LazyTB: Extended LTB-solutions with Multiple Interacting Fluids.
%                  Kimmo Kainulainen, Enrico Schiappacasse, Linda Tenhu, Olli Väisänen,
%                  XX.10.2026, arXiv: XXXX.XXXXX.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sol = NLTBFluid(input)



    yin = input.y(1);    % Starting point of evolution for solver parameters

    % Settings for pdepe.
    options = odeset('RelTol',input.Rtol,'AbsTol',input.Atol,'InitialStep',yin/1000);
   
    % Compute initial e-perturbation and create an interpolant. This is necessary
    % to compute the CDM-density perturbation inside the evolution equations.
    initial_e = initialize_runtime_variables(input);
    extra.initial_e = griddedInterpolant(input.tr,initial_e);

    % Run pdepe.
    ICon = @(r) InitialConditions(r,input);
    dfdx = @(r,y,u,DuDr) EvolutionEquations(r,y,u,DuDr,input,extra);
    sol  = pdepe(2,dfdx,ICon,@BoundaryConditions,input.tr,input.y,options);
    
end


%---------------------------------------------------------------------------------------------------
% FUNCTION defining the intial conditions
%---------------------------------------------------------------------------------------------------

function u0 = InitialConditions(r,input)
    
    ain = exp(input.y(1));                    % Initial scale factor.
    Nfluids = input.Nfluids;                   % Number of extra fluid components in addition to CDM
    
    
    % Initial density perturbation functions f and g  
    gri = zeros(1,Nfluids+1);
    fri = zeros(1,Nfluids+1);
    for ifluid = 1:Nfluids+1
        gri(ifluid) = input.gr_vec{ifluid}(r);
        fri(ifluid) = input.fr_vec{ifluid}(r);
    end

   
    % For the initial conditions we only need background EOS parameters,
    % and so I can set the Delta-vector here to 1, for all r.
    Deltas = ones(input.Nfluids+1,1);
                                    
    % Get EOS parameters
    [~,~,~,~,wbgV,wintV] = input.EOSPMS(ain,Deltas,[1,1]);

    %wbgV  = reshape(wbgV,[Nfluids+1,1]);
    %wintV = reshape(wintV,[Nfluids+1,1]);

    % Add CDM EOS parameters to the relevant vectors
    wbgV    = [0;wbgV];
    wint_DE = log(ain);
    wintV   = [wint_DE;wintV];

    % Calculate the mass fractions at the initial time. Note that this sum
    % includes CDM, vectors are length Nfluids+1!
    tOIin = input.Omega_vec'.*exp(-3.*wintV);          % Omega(a) for each species
    hbin  = sqrt(sum(tOIin));                          % hb == H/H0
    oIin  = tOIin/hbin.^2;                             % bar rho_i/bar rho; made a column vector!
    
    
    % Initial density contrast for dynamical (NON-CDM) fluids. Vector
    % size is (Nfluids,1).
    Delta_inV = (1 + fri(2:end)).';

    % Initials for e and dh. 
    RiV  = RiVec(oIin,wbgV);

    initial_e = - (ain*hbin)^2 * sum(oIin.*RiV.*gri.')/2;
    initial_dh = sum(oIin.*(1-RiV).*gri.')/2;

    % Non-trivial initial conditions for the velocity could be easily implemented 
    % by modifying this line.
    initial_v = zeros(input.Nfluids,1);

    % Initial conditions vector. Length 1+1+1+Nfluids+Nfluids.
    u0 = [ 0
           initial_e
           initial_dh
           log(Delta_inV) 
           initial_v ];


end

%------------------------------------------------------------------------
% Ri-integral, see Eq. (2.51) in [1]. -----------------------------------
%------------------------------------------------------------------------
function RiV = RiVec(Omg,w)

    RiV = Omg*0;
    
    DENOM = integral(@(u) integrand_RiVec(Omg,w,1,0,u) ,0,1);
    
    for i=1:length(Omg) 
        RiV(i) = integral(@(u) integrand_RiVec(Omg,w,i,1,u) ,0,1)/DENOM;
    end
    
end
function intgd = integrand_RiVec(Omg,w,i,k,u) 
    
    [U,OMG] = meshgrid(u,Omg);
    [~,  W] = meshgrid(u,w);
    
    DEN   = sum(OMG.*U.^(1-3*W),1).^1.5;
    
    intgd = u.^(3-k*(1+3*w(i)))./DEN;

end


%---------------------------------------------------------------------------------------------------
% FUNCTION defining the boundary conditions
%
% Direclet conditions 
%   etat --> 1;   e --> 0;    ht --> 1     DeW --> 1;   vp  --> 0
% vonNeumann conditions for derivatives cannot be used wVth LTB very easily
%---------------------------------------------------------------------------------------------------


function [pl,ql,pr,qr] = BoundaryConditions(~,~,~,ur,~)
    
    
    pl = zeros(length(ur),1);            % These are not relevant; because of the spherical symmetry
    ql =  ones(length(ur),1);            % the integrator ignores conditions set at r=0 boundary. It
                                         % adjusts them it internally.
    pr = ur;
    qr = zeros(length(ur),1);
    
    
end


%---------------------------------------------------------------------------------------------------
% FUNCTION defining the evolution equations for the metric
%---------------------------------------------------------------------------------------------------


function [c,f,s] = EvolutionEquations(r,y,u,DuDr,input,extra)

    Nfluids = input.Nfluids; % Number of non-CDM fluids.

    % Note: Suffixes below: V -> vector, r -> r-derivative
    
    % See [1] appendix A for definitions.
    detat  = u(1);                             % Perturbation in eta_t
    etatr  = DuDr(1);                          % r-derivative of detat
    e      = u(2);                             % Normalized curvature parameter
    dbt    = u(3);                             % Perturbation in b_t
    btr    = DuDr(3);                          % r-derivative of dbt  
    LDeV   = u(4:1:4+Nfluids-1);               % ln of density contrast for non-CDM fluids
    LDeVr  = DuDr(4:1:4+Nfluids-1);            % r-derivative of LDeV
    vpV    = u(4+Nfluids:1:4+2*Nfluids-1);     % velocities of non-CDM fluids
    vprV   = DuDr(4+Nfluids:1:4+2*Nfluids-1);  % r-derivative of vpV
    
    % Convert from perturbations to full quantities
    etat  = 1 + detat;                
    bt    = 1 + dbt;
    DeV   = exp(LDeV);            
    
    ab    = exp(y);            % Background scale factor
    
    % Find the cold dark matter density. 
    initial_e = extra.initial_e(r);
    initial_Delta_CDM = 1 + input.fr_vec{1}(r);    
    Delta_CDM = compute_Delta_CDM(initial_Delta_CDM,detat,etatr,initial_e,e,r);



    % NOTE: The fluid index should always be the FIRST dimension of Deltas!
    Deltas = [Delta_CDM, DeV.']; % Here fluid number is the last index.
    nd = ndims(Deltas);
    ord = 1:nd;
    Deltas = permute(Deltas,[nd,ord(1:end-1)]);


    % Define fluid equation of state parameters for each fluid componenet ..............................
    % NOTE: If you wish to optimize this program further, you should
    % consider writing your own version of NfluidEOSPMS and hardcoding the
    % fluids in your model.
    [cs2V,wV,wdotV,wrV,wbgV,wintV] = input.EOSPMS(ab,Deltas,[1,1]);


    % Out of the EOS parameters, we need to include CDM in wintV. A truncated 
    % version of wbgV is still necessary, so we save it for later use.
    wint_DM = log(ab);
    wintV   = [wint_DM; wintV];

    % The default arguments for the interaction terms also unclude the full 
    % EOS parameters w and cs2.
    wfull = [0;wV].';
    cs2full = [0;cs2V].';


    % Compute Omega(y) and scaled background Hubble rate.
    tOV = input.Omega_vec'.*exp(-3.*wintV);           % Omega(a) for each species
    hb  = sqrt(sum(tOV));                             % hb == H/H0
    oV  = tOV(2:end)/hb.^2;                           % bar rho_i/bar rho; made a column vector! NOTE: Only contains non-CDM-components.
    

    
    % Work out the various  coefficients in eos........................................................
    % See [1] appendix A.1. for the equations.
    
    % Lorentz-factors gamma_v^2 and gamma_cv^2
    gpV2  = 1./(1 - vpV.^2);
    gcvV2 = 1./(1 - cs2V.*vpV.^2);
    
    % Metric perturbation quantities. See appendix A.1 in [1].
    X   = sqrt(1 + 2*e*r^2);           
    etar = (etat + r*etatr)/X;                          
    chir = ab*hb*etar;
    chit = ab*hb*etat;
    
    % The tree sums involving EOS parameters in eq. (A.3).
    S1 = sum( oV.*(1+wV).*DeV.*vpV.*gpV2 );
    %S2 = sum( oV.*wbeffV );
    S2 = sum( oV.*wbgV );
    S3 = sum( oV.*(wV + vpV.^2).*gpV2.*DeV );
    
    % Compute b_r using the momentum constraint.
    dbr = dbt + (r*etat/(etat+r*etatr))*( btr - (3/2)*chir*S1 );         
    % The perturbation in the scaled local expansion rate.
    db  = (dbr + 2*dbt)/3; 
    
    % Convert perturbations to full quantities
    br  =  1 + dbr;   b  =  1 + db;
    

    % See Eq. (A.6) in [1].
    %nu = (1/3)*gcvV2.*( - vpV.^2.*(br - 3*cs2V*b) + 2*vpV*X/(r*chit) + vprV/chir );
    %nu = (1/3)*gcvV2.*( vpV.^2.*(br - 3*cs2V*b) + 2*vpV*X/(r*chit) + vprV/chir );

    %nu2 = 1 - b + (1/3)*gcvV2.*( vpV.^2.*(br - 3*cs2V*b) - 2*vpV*X/(r*chit) - vprV/chir );
    nu2 = -db + (1/3)*gcvV2.*( vpV.^2.*(br - 3*cs2V*b) - 2*vpV*X/(r*chit) - vprV/chir );
    %nu2 = 1 - gcvV2.*( b - vpV.^2.*br/3 + 2*vpV*X/(3*r*chit) + vprV/(3*chir)  );

    % The brackets in the v:derivative of Eq. (A.3) in [1].
    delA = - br + 3*cs2V*b + 2*cs2V.*vpV*X/(r*chit) - (1-cs2V).*vprV.*gpV2/chir;

    % See Eq. (A.5) in [1].
    B    = gcvV2.*(vpV.*wdotV + wrV./chir + cs2V.*LDeVr.*(1-vpV.^2)/chir);
    

    % ----------------------------------------------------------------------
    % Interaction terms ----------------------------------------------------
    % ----------------------------------------------------------------------

    vs = [0,vpV.']; 
    [Qt,arQr] = input.interactions(ab,vs,Deltas,wfull,cs2full);
        
    intPrefactor = gcvV2./ (DeV.*oV.*hb.^3);

    velInt = intPrefactor .* ( (1 + cs2V).*vpV.*Qt - (1 + cs2V.*vpV.^2).*arQr );
    RhoInt = intPrefactor .* ( (1 + vpV.^2).*Qt - 2*vpV.*arQr );


    % ----------------------------------------------------------------------
    % Right hand sides of the velocity and density equations ---------------
    % ----------------------------------------------------------------------

    % Velocity RHS
    velS = vpV.*(gcvV2./gpV2).*delA - (B + velInt)./(1+wV)./gpV2;

    % Density contrast RHS
    %RhoS = 3*( wbgV - wV - (1+wV).*(db + nu) ) - vpV.*LDeVr/chir + vpV.*B + RhoInt;

    RhoS = 3*( wbgV - wV + (1+wV).*nu2 ) - vpV.*LDeVr/chir + vpV.*B + RhoInt;



    % -------------------------------------------------------------------
    % Miscellaneous stuff -----------------------------------------------
    % -------------------------------------------------------------------

    % Freeze Delta for chosen fluids at a > input.freeze_scale_factor
    frz = (ab < input.freeze_scale_factor) .* input.freeze;


    % Check if we have a cosmological constant. In that case 1+vW=0 and as a result AV will be NaN's.
    % These must be converted to a zero to keep v_DE=0and rho_DE=const.
    IndexB = isnan(velS);  velS(IndexB) = 0; 

    
    % Finally add numerical diffusion regulators to damp perturbations at large r outside the void .......................................................
    % These depend on the application. Generally, large voids do not
    % necessarily need diffusion, but smaller ones do.
    
    Dd = input.DDiff(1/ab-1,r);   
    Dv = input.vDiff(1/ab-1,r);   
     

    % Eventual evolution equations using the above defined quantities..................................
    % See documentation of pdepe for the definitions of c,f,s.
    % --------------------------------------------------------------------------
    
    c = ones(3+2*Nfluids,1);    

    f = [0; 0; 0; frz.*Dd.*LDeVr; frz.*Dv.*vprV];    

    s = [  dbt*etat                                                     
           (3/2)*chit*X*S1/r
          -(3/2)*(bt*dbt - bt*S2 + S3) + e/chit^2 
           frz.*RhoS
           frz.*velS ];
    
    
end


%
% Calculate the initial value of e to an array. This is only done once and
% it is used to compute Delta_CDM inside the solver.
%
function initial_e = initialize_runtime_variables(input)
    initial_e = zeros(input.nr,1);
    for i = 1:input.nr
        u0 = InitialConditions(input.tr(i),input);
        ein = u0(2);
        initial_e(i) = ein;
    end
end

%
% Calculates the density contrast for pressureless matter
%
function Delta_CDM = compute_Delta_CDM(initial_Delta_CDM,detat,etatr,initial_e,e,r)
    
    etat = 1 + detat;
    Delta_CDM = initial_Delta_CDM./(etat.^2.*(etat + r.*etatr)).*sqrt((1 + 2*r.^2.*e)./(1 + 2*r.^2.*initial_e));
    
end


%---------------------------------------------------------------------------------------------------
% END CODE UNIT
%---------------------------------------------------------------------------------------------------
