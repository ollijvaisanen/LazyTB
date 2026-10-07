function [Delta_CDM,br_mom,br_ham,hb,b] = get_derived_quantities(solution,input)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    % Grid parameters: scaled radius and background scale factor.
    tr = solution.tr;
    dr = tr(2) - tr(1);

    sz = size(solution.detat);
    ny = sz(1);
    nr = sz(2);

    y = solution.y;
    ab = exp(y(1:ny));

    % Solution vectors
    bt   = 1 + solution.dbt;
    etat = 1 + solution.detat;
    e    = solution.e;


    % 4th order finite differences.
    etat_r = finite_diff_r_O4(etat,dr);
    e_r    = finite_diff_r_O4(e,dr);
    bt_r   = finite_diff_r_O4(bt,dr);



    % Collect fluid density contrasts and velocities
    fluid_names = input.fluid_names;
    Nfluids = input.Nfluids;
    Deltas = zeros(ny,input.nr,Nfluids+1);
    velocities = zeros(ny,input.nr,Nfluids+1);



    for i = 1:Nfluids
        fldname = sprintf("logDelta_%s",fluid_names{i});
        fldname2 = sprintf("v_%s",fluid_names{i});
        Deltas(:,:,i+1) = exp(solution.(fldname));
        velocities(:,:,i+1) = solution.(fldname2);
    end

    % Compute g_CDM(r) and its r-derivative.
    gr_CDM = input.gr_vec{1}(tr);
    gr_CDM_r = finite_diff_r_O4(gr_CDM,dr);


    Delta_CDM_init = 1 + input.gr_vec{1}(tr) + tr.*gr_CDM_r/3;
    e_init         = e(1,:);

    % Output this.
    Delta_CDM = get_Delta_CDM(tr,Delta_CDM_init,e_init,e,etat,etat_r);    
    
    % Permute Deltas to be compatible with the EOSPMS functions. 
    Deltas(:,:,1) = Delta_CDM;
    Deltas_EOSPMS = permute(Deltas,[3,1,2]);
    % Get EOS parameters
    [~,wV_Nfluids,~,~,~,wintV_Nfluids] = input.EOSPMS(ab,Deltas_EOSPMS,[ny,nr]);

    % Add CDM EOS parameterst to the arrays.
    wV_CDM    = zeros([1,ny,nr]);
    wintV_CDM = reshape(repmat(reshape(log(ab),[ny,1]),[1,nr]),[1,ny,nr]); % ffs


    wV    = [wV_CDM; wV_Nfluids];
    wintV = [wintV_CDM; wintV_Nfluids];

    % Cyclically permute the "fluid flavor" index to the last place.
    wV    = permute(wV,[2,3,1]);
    wintV = permute(wintV,[2,3,1]);

    % Stretch some arrays to make broadcasts more transparent.
    ab_rep = repmat(ab.',[1,nr,Nfluids+1]);
    tr_rep = repmat(tr,[ny,1]);

    % We need this later for Hamiltonian constraint
    X = sqrt(1 + 2*e.*tr_rep.^2);
    etar = (etat + tr_rep.*etat_r)./X;
  

    % Reshape Omega_vec ab for cleaner vectorized calculations, the
    % resulting shape is (ny,nr,Nfluids+1).
    Omega_vec = input.Omega_vec;
    Omega_vec_rep = permute(repmat(Omega_vec.',[1,ny,nr]),[2,3,1]);


    tOV = Omega_vec_rep.*exp(-3.*wintV);
    hb = sqrt(sum(tOV,3));
    oV = tOV./hb.^2;


    chit = ab_rep(:,:,1).*hb.*etat;  % In the paper notation, chi_\theta
    chir = ab_rep(:,:,1).*hb.*etar;  %     -||-          , chi_r
    gpV2 = 1./(1 - velocities.^2);


    sum1 = sum(oV.*(1+wV).*Deltas.*velocities.*gpV2,3);
    sum2 = sum(oV.*(1 + wV.*velocities.^2).*gpV2.*Deltas,3);
  
    br_mom = bt + (tr_rep.*etat./(etat+tr_rep.*etat_r)).*( bt_r - (3/2).*chir.*sum1 );         % Momentum constraint
    br_ham = -bt/2 + ((2*e + tr_rep.*e_r)./(chit.*chir.*X) + e./chit.^2)./bt + 3./(2*bt).*sum2;  % Hamiltonian constraint

    % Local Hubble rate
    b = (br_mom + 2*(1 + solution.dbt))/3;
 
    % hb is a background quantity and the same at all points.
    hb = hb(:,1);
end

function arr_r = finite_diff_r_O4(arr,dr)

    [ny,nr] = size(arr);
    slc = 3:nr+2;

    arr_pad = zeros(ny,nr+4);
    arr_pad(:,slc) = arr;
    arr_pad(:,2) = arr_pad(:,4);
    arr_pad(:,1) = arr_pad(:,5);
    arr_pad(:,nr+3) = arr_pad(:,nr+2);
    arr_pad(:,nr+4) = arr_pad(:,nr+2); 

    arr_r = (arr_pad(:,slc-2)/12 - 2*arr_pad(:,slc-1)/3 + 2*arr_pad(:,slc+1)/3 - arr_pad(:,slc+2)/12)/dr;
end

%
% Compute the density contrast for cold dark matter analytically.
%
function Delta_CDM = get_Delta_CDM(tr,Delta_CDM_init,e_init,e,etat,etat_r)

    % Reshape everything to correct shape for broadcasting
    sz = size(etat); na = sz(1); 

    tr = repmat(tr,[na,1]);
    Delta_CDM_init = repmat(Delta_CDM_init,[na,1]);
    e_init = repmat(e_init,[na,1]);
    
    Delta_CDM = Delta_CDM_init .* sqrt((1 + 2*tr.^2 .* e)./(1 + 2*tr.^2 .* e_init)) ./ etat.^2 ./ (etat + tr .* etat_r);
end