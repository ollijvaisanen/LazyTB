%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% A function for creating the initial conditions struct for the 
% light propagation solver.
%
%     [1]: LazyTB: Extended LTB-solutions with Multiple Interacting Fluids.
%                  Kimmo Kainulainen, Enrico Schiappacasse, Linda Tenhu, Olli Väisänen,
%                  XX.10.2026, arXiv: XXXX.XXXXX.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LP_input = create_light_prop_input(metric,LP_params)

    % Input struct for the light propagation function
    LP_input = struct();


    % ------------------------------------------------------------------------
    % ------ The background --------------------------------------------------
    % ------------------------------------------------------------------------

    a_grid = exp(metric.y);
    r_grid = metric.tr;


    % Create interpolants of the metric components
    method = "spline";         % As a_bg is not uniformly spaced, I'll use a spline. Alternatively I could make the interpolant in y = log(a_bg);
    extrapolation = "nearest"; % Here I'll give the fastest extrapolation method just to make the extrapolated values not NaN, and possible to set to zero by mutltiplication.
 
    Fdetat_interp = griddedInterpolant({a_grid,r_grid},metric.detat,method,extrapolation);
    Fe_interp     = griddedInterpolant({a_grid,r_grid},metric.e    ,method,extrapolation);

    
    rmax = max(r_grid);

    % Set detat and e to zero for r > trb.
    LP_input.Fdetat = @(a,r) (r < rmax) .* Fdetat_interp(a,r);
    LP_input.Fe     = @(a,r) (r < rmax) .* Fe_interp(a,r);


    % --------------------------------------------------------------------
    % Background scale factor and its time derivative --------------------
    % --------------------------------------------------------------------

    LP_input.adot    = metric.adot;
    LP_input.adotdot = metric.adotdot;


    % Perturbation in adot, adot = adot_bg*b.
    LP_input.b = griddedInterpolant({a_grid,r_grid},metric.b);

    % ------------------------------------------------------------------
    % Copy observer parameters from LP_params --------------------------
    % ------------------------------------------------------------------

    LP_input.ainit = LP_params.ainit;
    LP_input.r0    = LP_params.r0;
    LP_input.phi0  = LP_params.phi0;

    LP_input.freq0 = LP_params.freq0;
    
    LP_input.angles = LP_params.angles;
    LP_input.Nrays  = LP_params.Nrays;

    % --------------------------------------------------------------------
    % Copy other parameters ----------------------------------------------
    % --------------------------------------------------------------------

    LP_input.dr = LP_params.dr;
    LP_input.da = LP_params.da;

    % The end scale factor for light rays.
    LP_input.aend = LP_params.aend;

    % The end scale factor for the semi-analytic expansion evolution.
    LP_input.aend_semi_analytic = LP_params.factor_semi_analytic * LP_input.ainit;

    % The number of steps to save in the ODE45 solution.
    LP_input.Nsteps = LP_params.Nsteps;
    LP_input.Nsteps_semi_analytic = LP_params.Nsteps_semi_analytic;

    % Settings for the ode solver
    LP_input.settings = LP_params.settings;


    % -------------------------------------------------------------------
    % Compute derived quantities from the input parameters. -------------
    % -------------------------------------------------------------------
    
    LP_input = update_light_prop_input(LP_input);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function contains the physics part of the initial conditions. 
% Set the initial momentum vectors according to the given angles and
% metric, and get orthogonal frame vectors using Gram-Schmidt.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LP_input = update_light_prop_input(LP_input)
    
    r0 = LP_input.r0;
    dr = LP_input.dr;
    phi0 = LP_input.phi0;

    ainit = LP_input.ainit;


    % The metric components at the observer point.
    detat0 = LP_input.Fdetat(ainit,r0);
    e0     = LP_input.Fe(ainit,r0);

    % Differentiate detat wrt. (scaled) r at the observer
    ddetat0 = ( 1/12 * LP_input.Fdetat(ainit,r0 - 2*dr) - 2/3 * LP_input.Fdetat(ainit,r0 - dr) + 2/3 * LP_input.Fdetat(ainit,r0 + dr) - 1/12 * LP_input.Fdetat(ainit,r0 + 2*dr) )/dr;    

    Htt0 = ainit*r0 * (1 + detat0);
    Hrr0 = ainit * ( 1 + detat0 + r0 * ddetat0) ./ sqrt(1 + 2*e0*r0^2);


    % Observed angles for the light rays in radians. Angle = 0 points directly towards the origin.
    angles = LP_input.angles;

    LP_input.k_phi = zeros(LP_input.Nrays,1);

    for iray = 1:LP_input.Nrays

        k_r0   = LP_input.freq0 * Hrr0 * cos(angles(iray));
        k_phi0 = LP_input.freq0 * Htt0 * sin(angles(iray));



        % Set the initial state vector
        LP_input.yinit(iray,1) = r0;
        LP_input.yinit(iray,2) = phi0;
        LP_input.yinit(iray,3) = k_r0;

        % Set the (constant) angular momentum component
        LP_input.k_phi(iray) = k_phi0;

        % s2 frame vector using Gram-Schmidt. If there are numerical issues
        % change this to "if k_phi0 < some threshold".
        if k_phi0 ~= 0
            s2_normalization = 1/sqrt(Hrr0^2 + Htt0^2 * (k_r0 / k_phi0)^2);      
            s2r   =   s2_normalization;
            s2phi = - s2_normalization * k_r0 / k_phi0;
        else
            s2r = 0; 
            s2phi = 1/Htt0;
        end

        % The initial frame vectors
        LP_input.yinit(iray,4) = 0; % s2^0 = 0
        LP_input.yinit(iray,5) = s2r;
        LP_input.yinit(iray,6) = s2phi;

        % The expansion and shear parameters. 
        % These will be solved semi-analytically for the first couple
        % timesteps, so here they will be set to zero.

        expansion        = 0;
        sigma1           = 0;
        angular_distance = 0;

        LP_input.yinit(iray,7) = expansion;
        LP_input.yinit(iray,8) = sigma1;
        LP_input.yinit(iray,9) = angular_distance;


    end


    LP_input.angles = angles;

    % Reshape the initial state vector for the ODE-solver.
    LP_input.yinit = reshape(LP_input.yinit,LP_input.Nrays*9,1);
end