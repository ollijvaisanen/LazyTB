classdef NLTBSolver < handle
    %NLTB_SOLVER: This class handles the construction of the input struct for the 
    %             actual GR-solver, and packages the output in a readable form.
    % 
    %             This class contains no physics, and the solver function 
    %             NLTBFluid can be called directly if necessary. See
    %             the comments in NLTBFluid for the description of the input 
    %             structure.
    % 
    %             [1]: LazyTB: Extended LTB-solutions with Multiple Interacting Fluids.
    %                  Kimmo Kainulainen, Enrico Schiappacasse, Linda Tenhu, Olli Väisänen,
    %                  XX.10.2026, arXiv: XXXX.XXXXX.
    %
    %
    properties
        input
    end

    
    methods
        function obj = NLTBSolver(ka)
            %NLTB_SOLVER Constructor
            % ka = Inpt keyword arguments as NLTBSolver(h=0.76,...)
            arguments

                ka.h = 0.71;              % Hubble rate H0 in 100 Mpc/km/s.

                ka.trgd = 8 * 100 * 0.71 * 3.336e-4; % Default grid boundary at 800 Mpc/h, in units if c/H0.

                ka.zin  = 1000;          % Initial solver redshift
                ka.zout = 0;             % Final solver redshift

                % Default parameters for the reference CDM-fluid.
                ka.Omega_CDM = 1.0;                            % Background density fraction of cold dark matter
                ka.gr_CDM    = @(r) zeros(size(r));            % Initial CDM profile function g(r)   
                ka.fr_CDM    = @(r) zeros(size(r));            % f(r) = (1 + r/3 * d/dr)g(r).


                % Default fluid parameters contain no additional fluids.
                ka.Omegas       = [];
                ka.grs          = {};
                ka.frs          = {};
                ka.fluid_names  = {};

                % If no other fluids are declared, this should just return empty arrays.
                ka.EOSPMS = @(a,Deltas,sz) deal([],[],[],[],[],[]);

                ka.nr = 150;     % Number of grid points for radius r;
                ka.ny = 50;      % Number of grid points for log of scale factor y = log(a/a0);

                % Solver tolerances. This code uses the matlab solver
                % pdepe, see its documentation for details on solver
                % parameters.
                ka.Rtol = 1e-13; % Relative solver tolerance.
                ka.Atol = 1e-13; % Absolute solver tolerance.

                % Default numerical diffusion set to zero.
                ka.DDiff = @(z,tr) zeros(size(tr));
                ka.vDiff = @(z,tr) zeros(size(tr));

                % Default interaction terms set to zero.
                ka.interactions = @(a,v,Delta,w,cs2) deal(0,0);
                
                % Optionally allow the user to freeze one or more fluids
                % for a > freeze_scale_factor. Off by default.
                ka.freeze = 1;
                ka.freeze_scale_factor = inf;

                % If you wish to construct the input struct directly, any
                % parameters set to this argument will override their
                % equivalents.
                ka.input = NaN;

            end
            
            % Set the solver parameters. The below code overrides the above default settings.

            % The solver functions take the parameters inside a struct.
            obj.input.h     = ka.h;
            obj.input.trgd  = ka.trgd;
            obj.input.zin   = ka.zin;
            obj.input.zout  = ka.zout;
            obj.input.nr    = ka.nr;
            obj.input.ny    = ka.ny;
            obj.input.Rtol  = ka.Rtol;
            obj.input.Atol  = ka.Atol;
            obj.input.DDiff = ka.DDiff;
            obj.input.vDiff = ka.vDiff;

            % Compute the calculation grids.
            yin            = -log(1+obj.input.zin);                   % Initial value of the integration variable y = log a
            yout           = -log(1+obj.input.zout);                  % Initial value of the integration variable y = log a
            obj.input.tr   = linspace(0,ka.trgd,ka.nr);               % grid in scaled r-variable
            obj.input.y    = linspace(yin,yout,ka.ny);                % log of background scale factor


            % Setup the fluid parameters ---------------------------------

            obj.input.Nfluids = length(ka.Omegas);
            
            obj.input.Omega_vec = zeros(1,obj.input.Nfluids);              % Omega_vec also contains CDM as its first component
            obj.input.Omega_vec(1) = ka.Omega_CDM;                         % We need to make sure Omega_vec is a (1,Nfl+1)-array. 
            obj.input.Omega_vec(2:obj.input.Nfluids+1) = ka.Omegas(:);     

            obj.input.fluid_names = ka.fluid_names;

            obj.input.gr_vec = {ka.gr_CDM};
            obj.input.fr_vec = {ka.fr_CDM};
            for i = 1:obj.input.Nfluids
                obj.input.gr_vec{end+1} = ka.grs{i};
                obj.input.fr_vec{end+1} = ka.frs{i};
            end

            obj.input.EOSPMS = ka.EOSPMS;

            % Set the interactions between the fluids --------------------
            obj.input.interactions = ka.interactions;

            % Allow the option to freeze one or more fluid equations -----
            obj.input.freeze = ka.freeze;
            obj.input.freeze_scale_factor = ka.freeze_scale_factor;

            % If an input struct is passed directly, the fields in this
            % struct will override those set above.
            if isstruct(ka.input)
                input_fldnames = fieldnames(ka.input);
                for i = 1:length(input_fldnames)
                    obj.input.(input_fldnames{i}) = ka.input.(input_fldnames{i});
                end
            end


            
        end

  

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Running the simulation
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
        function output = run(obj,ka)

            arguments
                obj
                ka.add_input = true;   % Choose if you want to append input struct to output.
                ka.add_derived = true; % Choose if you want to compute additional derived quantities. See below.
            end

            % Run the solver
            tic
            solN    = NLTBFluid(obj.input);
            toc

            % --------------------------------------------------------------------
            % Repackage the solution into a more readable form -------------------
            % --------------------------------------------------------------------
    
            % NOTE: See appendix A.1 in [1] for the explanations of variables.    
            output.solution.detat = solN(:,:,1);
            output.solution.e = solN(:,:,2);
            output.solution.dbt = solN(:,:,3);
    
            for i = 1:obj.input.Nfluids

                namestr_DeV = sprintf("logDelta_%s",obj.input.fluid_names{i});
                namestr_vpV = sprintf("v_%s",obj.input.fluid_names{i});

                output.solution.(namestr_DeV) = solN(:,:,3+i);
                output.solution.(namestr_vpV) = solN(:,:,3+obj.input.Nfluids+i);       
             end
            
            % Add the grids to the solution struct.
            output.solution.tr = obj.input.tr;
            output.solution.y  = obj.input.y;

            if ka.add_input
                % Append the full input struct to the output for chekups.
                output.input = obj.input;
            end

            if ka.add_derived
                % Compute derived quantities and add to the output. See below.
                output.derived = get_derived_quantities(obj,output.solution);
            end
        end    

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Compute derived quantities from the solution
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function out = get_derived_quantities(obj,sol)
            
            [Delta_CDM,br_mom,br_ham,hb,b] = get_derived_quantities(sol,obj.input);

            out.br_mom = br_mom;              % hr via momentum constraint
            out.br_ham = br_ham;              % hr via hamiltonian constraint
            out.Delta_CDM = Delta_CDM;        % Density contrast of dust
            out.hb = hb;                      % H/H0, H = hubble rate
            out.b = b;                        % Local Hubble rate, b = (b_r + 2*b_theta)/3. See eq. (A.1) in [1].
        end
    end
end

