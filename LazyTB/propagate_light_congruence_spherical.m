%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% A function for calculating the trajectories of light rays in a
% spherically symmetric spacetime. Also solves the Sachs equations for the
% expansion and shear parameters.
%
% The input struct has the following fields. See appendix A.1 in [1] for 
% definitions:
%
%       Fdetat: Function handle returning detat(a,r).
%           Fe: Function handle returning e(a,r).
%         adot: Function handle returning adot(a) in units of H0.
%      adotdot: Function handle returning adotdot(a) in units of H0^2.
%            b: Function handle returning b(a,r).
%        ainit: Scale factor at the observer timeslice.
%         aend: Scale factor at the end of run. Note that aend < ainit.
%           r0: Observer radius.
%         phi0: Observer angle.
%        freq0: 0-coordinate of outgoing light rays at the start. In a run 
%               with no numerical error this has no effect on output.
%       angles: (1,Nrays) array. Angles of outgoing rays, with 0 pointing
%               directly at the origin.
%           dr: r-step used internally in finite differences.
%           da: a-step used internally in finite differences.
%        Nrays: Number of outgoing rays.
%   aend_semi_analytic: End scale factor of the semi-analytic evolution for
%               setting the initial conditions to expansion scalar and DA.
%               The default value is 0.99999*ainit, a too large value will
%               usually result in a constant shift of DA.
%       Nsteps: Number of output steps for ode45.
%  Nsteps_semi_analytic: Number of semi-analytic timesteps for expansion 
%               and DA. 20 by default.
%     settings: Standard settings object for ode45.
%        k_phi: Initial angular momentum for outgoing rays. This was set by
%               the solver object based on the angles.
%        yinit: ln(ainit).
%
%  Note that many of the above are related to each other and they are set
%  automatically by the LightProbSolver.
%
%     [1]: LazyTB: Extended LTB-solutions with Multiple Interacting Fluids.
%                  Kimmo Kainulainen, Enrico Schiappacasse, Linda Tenhu, Olli Väisänen,
%                  XX.10.2026, arXiv: XXXX.XXXXX.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = propagate_light_congruence_spherical(input)



    % --------------------------------------------------------------------
    % Evolve expansion and shear analytically for a while ----------------
    % --------------------------------------------------------------------
    
    uspan_semi_analytic = linspace(0,input.ainit-input.aend_semi_analytic,input.Nsteps_semi_analytic);
   

    [u_semi_analytic,y_semi_analytic] = ode45(@(u,y) evolution_equations(u,y,input,true),uspan_semi_analytic,input.yinit,input.settings);
    

    y_semi_analytic = reshape(y_semi_analytic,input.Nsteps_semi_analytic,input.Nrays,9);
    
    y_semi_analytic = fill_expansion_angular_distance(u_semi_analytic,y_semi_analytic,input);
    
    % Get the initial state vector for the full numerical evolution.
    new_yinit = reshape(y_semi_analytic(end,:,:),input.Nrays,9);

    
    % --------------------------------------------------------------------
    % Full numerical evolution -------------------------------------------
    % --------------------------------------------------------------------
    
    
    uspan = linspace(0,input.aend_semi_analytic-input.aend,input.Nsteps);
    
    [u_numerical,y_numerical] = ode45(@(u,y) evolution_equations(u,y,input,false),uspan,new_yinit,input.settings);

    y_numerical = reshape(y_numerical,input.Nsteps,input.Nrays,9);
    
    y = [y_semi_analytic ; y_numerical(2:end,:,:)];
    
    % Geodesic equation variables. Note that kphi = constant.
    out.a = [input.ainit - u_semi_analytic; input.aend_semi_analytic - u_numerical(2:end)];
    out.rrays   = y(:,:,1);
    out.phirays = y(:,:,2);
    out.krrays  = y(:,:,3);
    
    % The congruence equation variables.
    out.s2        = y(:,:,4:6);
    out.expansion = y(:,:,7);
    out.sigma1    = y(:,:,8);
    out.angular_distance = y(:,:,9);
    
    out.k_phi = input.k_phi;
    
    
    [zrays,k0] = calculate_redshift(out,input); 
    
    out.zrays = zrays;
    out.k0    = k0;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The right hand side of the geodesic equations + sachs equations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydu = evolution_equations(u,y,input,semi_analytic)

    

    % Reshape y to separate different light rays.
    y = reshape(y,input.Nrays,9);
    
    % The background scale factor
    if semi_analytic
        a = input.ainit - u;
    else
        a = input.aend_semi_analytic - u;
    end
    a_arr = a*ones(input.Nrays,1); % For convenience.
    
    % Time derivatives of the BACKGROUND scale factor.
    adot = input.adot(a);
    adotdot = input.adotdot(a);
    

    
    % Unpack the state vector for readability. 
    rrays = abs(y(:,1)); % Take the absolute value of radius for the interpolation. THIS MIGHT CAUSE ISSUES NEAR THE ORIGIN!
    % phirays = y(:,2); % Due to spherical symmetry, phi is never needed explicitly
    krrays = y(:,3);
    
    %fprintf("u = %d, r = %d\n",u,rrays(1))
    
    % The frame vector s2. s1 will be set analytically after the metric
    % components have been calculated.
    s20 = y(:,4);
    s2r = y(:,5);
    s2p = y(:,6);
    
    
    % Expansion and shear tensor. Note that sigma2 = 0.
    expansion        = y(:,7);
    sigma1           = y(:,8);
    angular_distance = y(:,9);
    
    
    % The phi-component of momentum is a constant of motion due to
    % spherical symmetry.
    k_phi = input.k_phi;
    
    
    % NOTE: I could optimize this by tabulating values of the interpolants.
    % I'll do that once the program works.
    dr = input.dr;
    da = input.da;
    
    % Tabulate the detat- FD points for efficient vectorization.
    rpts = zeros(input.Nrays,5,5);
    apts = zeros(input.Nrays,5,5);
    for i = -2:2
        for j = -2:2
            rpts(:,3+i,3+j) = abs(rrays + j*dr);
            apts(:,3+i,3+j) = a_arr + i*da;
        end
    end
    reshape(rpts,25*input.Nrays,1);
    reshape(apts,25*input.Nrays,1);
    
    detat_vals = reshape(input.Fdetat(apts,rpts),input.Nrays,5,5);
    
    % I'll use the same tabulated evaluation points. However, I don't need
    % the "mixed points" for e, so I'll have to take some slices.
    rpts = reshape(rpts,input.Nrays,5,5);
    apts = reshape(apts,input.Nrays,5,5);
    
    % Tabulate e-values for different a.
    e_vals = reshape(input.Fe(reshape(apts(:,:,3),input.Nrays*5,1),reshape(rpts(:,:,3),input.Nrays*5,1)),input.Nrays,5);
    
    % Tabulate the e-values for different r.
    e_vals_r = reshape(input.Fe(reshape(apts(:,3,:),input.Nrays*5,1),reshape(rpts(:,3,:),input.Nrays*5,1)),input.Nrays,5);
    
    
    % The metric perturbations
    detat = detat_vals(:,3,3);
    e     = e_vals(:,3);

    
    % Calculate the r-derivatives of detat, e numerically to 4th order.
    % NOTE: This version will fail too close to the origin, though the
    % resulting dipoles still seem very close to the cartesian version.
    
      
        
    % d(detat)/dr
    ddetat_dr = ( 1/12 * detat_vals(:,3,1) - 2/3 * detat_vals(:,3,2) + 2/3 * detat_vals(:,3,4) - 1/12 * detat_vals(:,3,5) )/dr;    
   
    % d^2(detat)/dr^2
    d2detat_dr2 = ( -1/12 * detat_vals(:,3,1) + 4/3 * detat_vals(:,3,2) - 5/2 * detat_vals(:,3,3) + 4/3 * detat_vals(:,3,4) - 1/12 * detat_vals(:,3,5) )/dr^2;
   
    % TIME DERIVATIVES: Note: If necessary, I'll substitute the evolution
    % equations for some of these.
    
    % d(detat)/da
    ddetat_dt  =  adot * ( 1/12 * detat_vals(:,1,3) - 2/3 * detat_vals(:,2,3) + 2/3 * detat_vals(:,4,3) - 1/12 * detat_vals(:,5,3) )/da;
    
    % d(detat)/dt^2 
    d2detat_da2 = ( -1/12 * detat_vals(:,1,3)                    + 4/3 * detat_vals(:,2,3)              - 5/2 * detat_vals(:,3,3)         + 4/3 * detat_vals(:,4,3)              - 1/12 * detat_vals(:,5,3))/da^2;
    
    d2detat_dt2 = adotdot/adot * ddetat_dt + adot^2 * d2detat_da2;
    
    % d^2(detat)/dtdr
    
    d2detat_dtdr = adot * (      detat_vals(:,1,1)  -  8 * detat_vals(:,1,2)  +  8 * detat_vals(:,1,4) -     detat_vals(:,1,5)   ...
                            -8 * detat_vals(:,2,1)  + 64 * detat_vals(:,2,2)  - 64 * detat_vals(:,2,4) + 8 * detat_vals(:,2,5)   ...
                            +8 * detat_vals(:,4,1)  - 64 * detat_vals(:,4,2)  + 64 * detat_vals(:,4,4) - 8 * detat_vals(:,4,5)   ...
                            -    detat_vals(:,5,1)  +  8 * detat_vals(:,5,2)  -  8 * detat_vals(:,5,4) +     detat_vals(:,5,5)) / (144*da*dr);
    
    
    
    
    d3detat_da2dr = ( -     detat_vals(:,1,1) +   8 * detat_vals(:,1,2) -   8 * detat_vals(:,1,4) +      detat_vals(:,1,5) ...
                      +16 * detat_vals(:,2,1) - 128 * detat_vals(:,2,2) + 128 * detat_vals(:,2,4) - 16 * detat_vals(:,2,5) ...
                      -30 * detat_vals(:,3,1) + 240 * detat_vals(:,3,2) - 240 * detat_vals(:,3,4) + 30 * detat_vals(:,3,5) ...    
                      +16 * detat_vals(:,4,1) - 128 * detat_vals(:,4,2) + 128 * detat_vals(:,4,4) - 16 * detat_vals(:,4,5) ...
                      -     detat_vals(:,5,1) +   8 * detat_vals(:,5,2) -   8 * detat_vals(:,5,4) +      detat_vals(:,5,5) )/(144 * dr * da^2);
                                     
    d3detat_dt2dr = adotdot/adot * d2detat_dtdr + adot^2 * d3detat_da2dr;
    
    
    
    de_dt = adot * ( 1/12 * e_vals(:,1) - 2/3 * e_vals(:,2) + 2/3 * e_vals(:,4) - 1/12 * e_vals(:,5) )/da;
    
    d2e_da2 = ( -1/12 * e_vals(:,1) + 4/3 * e_vals(:,2) - 5/2 * e_vals(:,3) + 4/3 * e_vals(:,4) - 1/12 * e_vals(:,5))/da^2;
     
    d2e_dt2 = adotdot/adot * de_dt + adot^2 * d2e_da2;
    
    % de/dr
    de_dr           = ( 1/12 * e_vals_r(:,1) - 2/3 * e_vals_r(:,2) + 2/3 * e_vals_r(:,4) - 1/12 * e_vals_r(:,5) )/dr;
    
 
    % Calculate the full metric components and their derivatives.
      
    % NOTE: a_theta in [1] is same as F. Here we switched notation to keep 
    % the code at least somewhat clean.
    F        = a * ( 1 + detat ); 
    dF_dr    = a * ddetat_dr;
    d2F_dr2  = a * d2detat_dr2;
    dF_dt    = adot/a * F + a * ddetat_dt;
    d2F_dt2  = (adotdot/a - (adot/a)^2) * F + adot/a * dF_dt + adot * ddetat_dt + a * d2detat_dt2;
    d2F_dtdr = adot/a * dF_dr + a * d2detat_dtdr;
    
    %Htt = rrays .* F;
    %dHtt_dr = F + rrays .* dF_dr;
    
    Hrr       = (F + rrays .* dF_dr) ./ sqrt(1 + 2*e.*rrays.^2);
    dHrr_dr = a * ( 2*ddetat_dr + rrays .* d2detat_dr2 )./sqrt(1 + 2*e.*rrays.^2) - a/(2) * ( 1 + detat + rrays .* ddetat_dr) .* ( 4*rrays.*e + 2*rrays.^2 .* de_dr ) ./ ( 1 + 2 * rrays.^2 .* e ).^(3/2);   
  
    dHrr_dt   = (1 + 2*e.*rrays.^2).^(-3/2) .* ( (adot * ( 1 + detat + ddetat_dr .* rrays) + a*ddetat_dt + a*rrays .* d2detat_dtdr ) .* (1 + 2*e.*rrays.^2) - a * de_dt .* rrays.^2 .* ( 1 + detat + ddetat_dr .* rrays) );  
    d2Hrr_dt2 = (1 + 2*e.*rrays.^2).^(-5/2) .* ( 3*a*de_dt.^2 .* rrays.^4 .* (1 + detat + ddetat_dr .* rrays) - ( 2*a*de_dt .* rrays.^2 .* (ddetat_dt + d2detat_dtdr .* rrays) + 2*adot*de_dt .* rrays.^2 .* (1 + detat + ddetat_dr .* rrays) - a*d2e_dt2 .* rrays.^2 .* (1 + detat + ddetat_dr .* rrays) ) .* (1 + 2*e.*rrays.^2) ...
                                                       + ( 2*adot*(ddetat_dt + d2detat_dtdr .* rrays) + a*(d2detat_dt2 + d3detat_dt2dr .* rrays) + adotdot*(1 + detat + ddetat_dr .* rrays) ) .* (1 + 2*e.*rrays.^2).^2 );                 
    

    % The zero-component of the photon momentum (upper index)
    k0 = sqrt( krrays.^2 ./ Hrr.^2 + k_phi.^2 ./ (rrays .* F).^2 );
    

    % Set the first frame vector.
    s1t = - 1./(rrays .* F);
    
    dydu = zeros(input.Nrays,9);
    
    % The Jacobian factor in front of every evolution equation. Note that
    % df/da = jac df/du
    jac = (adot .* k0).^(-1);
    
    % --------------------------------------------------------------------
    % The geodesic equations --------------------------------------------
    

    % dr/da
    dydu(:,1) = - jac .* krrays ./ Hrr.^2;
    % dphi/da
    dydu(:,2) = - jac .* k_phi ./ (rrays .* F).^2;
    % dkr/da
    dydu(:,3) = - jac .* ( krrays.^2 ./ Hrr.^3 .* dHrr_dr  +  k_phi.^2 .* ( 1./(rrays.^3 .* F.^2) + dF_dr ./ (rrays.^2 .* F.^3) ) );

    % The parallel transport equations for s2. s1 is handled analytically
    
    % ds2^0/da
    dydu(:,4) = - jac .* ( - krrays .* s2r .* dHrr_dt ./ Hrr - k_phi .* s2p .* dF_dt ./ F );
    % ds2^r/da
    dydu(:,5) = - jac ./ Hrr.^2 .* ( - Hrr .* dHrr_dt .* k0 .* s2r - krrays ./ Hrr .* (dHrr_dt .* s20 + dHrr_dr .* s2r) + k_phi .* s2p .* (dF_dr ./ F + rrays.^(-1)) );
    % ds2^p/da
    dydu(:,6) = - jac .* ( - ( krrays .* dF_dr .* s2p ./ Hrr.^2 + dF_dt .* k0 .* s2p )./F - k_phi .* s2r ./ (F.^2 .* rrays.^3) - k_phi .* (dF_dt .* s20 + dF_dr .* s2r)./(F.^3 .* rrays.^2) - krrays .* s2p ./ (Hrr.^2 .* rrays) ); 
    
    if ~semi_analytic
        
        % Evolution equations for the expansion and shear.
        
        % Contraction Ricci_{mu nu} k^mu k^nu.
        Ricci_kk =   k0.^2        .* ( - 2 * d2F_dt2 ./ F - d2Hrr_dt2 ./ Hrr ) ...
                   + k_phi.^2     .* ( ( 1./F.^4 - 1./(F.^2 .* Hrr.^2) )./rrays.^4  + ( - 4 * dF_dr ./(F.^3 .* Hrr.^2) + dHrr_dr./(F.^2 .* Hrr.^3) )./rrays.^3   ...
                                     + ( dF_dt.^2 ./ F.^4 + dF_dt .* dHrr_dt ./ (F.^3 .* Hrr) + d2F_dt2 ./ F.^3 - dF_dr.^2 ./ (F.^4 .* Hrr.^2) + dF_dr .* dHrr_dr ./ (F.^3 .* Hrr.^3) - d2F_dr2 ./ (F.^3 .* Hrr.^2) )./rrays.^2 ) ...
                   + krrays .* k0 .* ( 4 * dHrr_dt .* dF_dr ./ (F .* Hrr.^3) - 4 * d2F_dtdr ./ (F .* Hrr.^2) + (-4 * dF_dt ./ (F .* Hrr.^2) + 4*dHrr_dt ./ Hrr.^3 )./rrays ) ...
                   + krrays.^2    .* ( 2 * dF_dt .* dHrr_dt ./ (F .* Hrr.^3) + d2Hrr_dt2 ./ Hrr.^3 + 2 * dF_dr .* dHrr_dr ./ (F .* Hrr.^5) - 2 * d2F_dr2 ./ (F.* Hrr.^4) + (-4*dF_dr./F + 2*dHrr_dr./Hrr)./(Hrr.^4 .* rrays) );
        
                     
        riemann_ricci_diff =   k0.^2     .* d2Hrr_dt2 ./ Hrr / 2 ...
                             - krrays.^2 .* d2Hrr_dt2 ./ Hrr.^3 / 2 ...
                             + k_phi.^2  .* ( (1./(2*F.^4) - 1./(2*F.^2 .* Hrr.^2))./rrays.^4 - dHrr_dr./(2 * F.^2 .* Hrr.^3 .* rrays.^3) ...
                                            + ( dF_dt.^2 ./ (2*F.^4) - dF_dt .* dHrr_dt ./ (2 * F.^3 .* Hrr) - d2F_dt2 ./ (2*F.^3) - dF_dr.^2 ./ (2 * F.^4 .* Hrr.^2) - dF_dr .* dHrr_dr ./ (2 * F.^3 .* Hrr.^3) + d2F_dr2 ./ (2 * F.^3 .* Hrr.^2) )./rrays.^2 );
        
    
        
        % d(expansion)/da
        dydu(:,7) = - jac .* ( - 1/2 * expansion.^2 - 2 * sigma1.^2 - Ricci_kk );

        % d(sigma1)/da
        %dydu(:,8) = - jac .* ( - expansion .* sigma1 - ( Riemann_sskk - Ricci_kk / 2 ) );
        dydu(:,8) = - jac .* ( - expansion .* sigma1 - riemann_ricci_diff );

        % d(angular_distance)/da
        dydu(:,9) = - jac .* expansion .* angular_distance / 2;
        

    end
    
    %dydu(:,9)
    
    dydu = reshape(dydu,input.Nrays*9,1);
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Fills the expansion, shear and angular distance of the first few
% timesteps under the assumption of a local FRW-universe.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function yout = fill_expansion_angular_distance(u,y,input)
    
    yout = y;

    a = input.ainit - u;
    rrays = y(:,:,1);
    sz = size(rrays); Nrays = sz(2);

    a_rep = repmat(a,[1,Nrays]);

    % Interpolate b for every ray at every timestep.
    brays = input.b(repmat(a,[1,Nrays]),rrays);

    % NOTICE! I set adot using an inhomogeneous value for H0 at the origin
    % points of the rays. This ignores the (slight) evolution in the
    % deviation of the adot from frw due to the bubble. 
    %adot  = repmat(input.adot(a),[1,Nrays]) .* hrays(1,:)/3;
    adot  = repmat(input.adot(a),[1,Nrays]) .* brays/3;
    


    D_phys_integrand = 1./(a_rep .* input.adot(a));
    %D_phys_integrand = 1./(a_rep .* adot);
    D_phys = cumtrapz(flip(a),flip(D_phys_integrand,1),1);
    D_phys = flip(D_phys(end) - D_phys,1);
    
    angular_distance = a_rep .* D_phys;

    %angular_distance
    %error
    
    k0_frw = input.freq0 ./ a;
    
    %expansion = 2 * k0_frw .* ( adot./a_rep - 1./angular_distance );
    expansion = 2 * k0_frw .* ( input.adot(a)./a_rep - 1./angular_distance );
    
    yout(:,:,7) = expansion;
    yout(:,:,9) = angular_distance;

    %yout(:,:,7) = repmat(reshape(expansion,Na,1),1,input.Nrays);
    %yout(:,:,9) = repmat(reshape(angular_distance,Na,1),1,input.Nrays);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  A function for calculating the redshift along a set of light
%  trajectories. For now, I'll assume that both the source and the observer
%  are stationary.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [zrays,k0_rays] = calculate_redshift(solution,input)
    
    Nsteps = length(solution.a);
    
    % For convenience.
    arays = repmat(solution.a,1,input.Nrays); 
    k_phi_rays = repmat(input.k_phi.',Nsteps,1);
    
    rrays = solution.rrays;
    
    % The metric perturbations
    detat = input.Fdetat(arays,rrays);
    e = input.Fe(arays,rrays);
    

    
    % Compute r-derivative of detat using finite differences.
    dr = input.dr;
    ddetat  = ( 1/12 * input.Fdetat(arays,rrays - 2*dr) - 2/3 * input.Fdetat(arays,rrays - dr) + 2/3 * input.Fdetat(arays,rrays + dr) - 1/12 * input.Fdetat(arays,rrays + 2*dr) )/dr;    
    
    % The metric components along the rays.
    Htt  = arays .* rrays .* ( 1 + detat ); 
    Hrr = arays .* ( ( 1 + detat) + rrays .* ddetat ) ./ sqrt(1 + 2*e.*rrays.^2);
    
    
    % Calculate k0 along the trajectories.
    k0_rays = sqrt(solution.krrays.^2 ./ Hrr.^2 + k_phi_rays.^2 ./ Htt.^2);
    
    % Get k0 measured by the observer and compute the observed redshifts.
    k0_obs = k0_rays(1,:);    
    zrays = k0_rays ./ repmat(k0_obs,Nsteps,1) - 1;
end