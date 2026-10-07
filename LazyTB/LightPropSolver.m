classdef LightPropSolver
    %UNTITLED Solver class for light propagation
    %   Detailed explanation goes here

    properties
        LP_input   % Iput struct to be given to the solver
        LP_params  % Parameter struct for constructing the input
        %NLTB_solution   % Solution and input structs from NLTB-solvers.

        metric % Necessary parameters from the metric

        LP_solution % Solution struct. This will be NaN until you run the solver for the first time.
    end

    methods
        function obj = LightPropSolver(metric,kw)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here;
            arguments
                %NLTB_solution 

                metric

                kw.ainit = 1.0;
                kw.r0    = 1e-5;
                kw.phi0  = 0.0;
    
                kw.freq0 = 1.0; % Initial dimensionless frequency. 
                                % Note that redshifts and distances are always independent of freq0,
                                % and this variable simply scales the internal numbers.

                kw.angles = [pi];

                kw.aend = 1/(1 + 10);

                % If no value is given for dr, use the solution grid spacing, see below.
                kw.dr = NaN; 
                kw.da = 1e-3;

                kw.Nsteps = 1000;
                kw.Nsteps_semi_analytic = 20;
                kw.factor_semi_analytic = 0.999999;
                
                % If no value is give, use default settings. See below.
                kw.ode_settings = NaN;
            end
            
            %obj.NLTB_solution = NLTB_solution;
            

            % -----------------------------------------------------------
            % Create light propagation parameters struct ----------------
            % -----------------------------------------------------------

            % Initial and final scale factors. Note that the solver integrates
            % backwards in time from the observer.
            LP_params.ainit = kw.ainit;
            LP_params.aend = kw.aend;

            % Initial radius and angle. The polar angle is assumed to be pi/2.
            LP_params.r0    = kw.r0;
            LP_params.phi0  = kw.phi0;

            LP_params.freq0 = kw.freq0;

            LP_params.angles = kw.angles;
            LP_params.Nrays  = length(kw.angles);

            if isnan(kw.ode_settings)
                % Default settings for the ode solver.
                LP_params.settings = odeset("RelTol",1e-10,"AbsTol",1e-10);
            else
                LP_params.settings = kw.ode_settings;
            end

            % Set the spacing for the finite differences used to differentiate the met
            LP_params.da = kw.da;
            if isnan(kw.dr)
                % If no value is given for dr, use the solution grid spacing.
                LP_params.dr = (NLTB_solution.tr(2) - NLTB_solution.tr(1));
            else
                LP_params.dr = kw.dr;
            end

            % Set the number of timesteps in the output.
            LP_params.Nsteps = kw.Nsteps;

            % In order to set the initial conditions for the angular distance,
            % part of the equations is computed analytically in the local FLRW
            % metric until scale factor factor_semi_analytic * ainit.
            LP_params.Nsteps_semi_analytic = kw.Nsteps_semi_analytic;
            LP_params.factor_semi_analytic = kw.factor_semi_analytic;
            
            % Creates an input struct out of the parameters. 
            obj.LP_input = create_light_prop_input(metric,LP_params);

            obj.LP_solution = NaN;


        end

        function output = run(obj)
            out = propagate_light_congruence_spherical(obj.LP_input);
    
            % Add angle values to the output struct for completeness.
            out.angles = obj.LP_input.angles;

            output = out;
        end

        % Computes the differences in output variables to a background
        % FLRW-model.
        function diff = get_output_differences(obj,sol)


            Nrays = length(obj.LP_input.angles);

            abg = sol.a;
            ainit = sol.a(1);

            zbg = 1./(abg/ainit) - 1;
            zbg_rep = repmat(zbg,1,Nrays);
            

            % Relative differences in redshift and observed temperature
            % from a FLRW-model.
            diff.zrays_relative = (sol.zrays - zbg_rep)./zbg_rep;
            diff.temperature_relative = (1 + zbg_rep)./(1 + sol.zrays) - 1;

            % Integrate background physical distance.
            D_phys_integrand = 1./(abg .* obj.LP_input.adot(abg));
            D_phys = cumtrapz(flip(abg),flip(D_phys_integrand));
            D_phys = flip(D_phys(end) - D_phys); 

            % Background angular diamtere distance.
            DA_bg = abg/ainit .* D_phys;
            DA_bg_rep = repmat(DA_bg,1,Nrays);

            diff.angular_distance_relative = (sol.angular_distance - DA_bg_rep) ./ DA_bg_rep;
            diff.angular_distance_bg = DA_bg_rep;

        end
    end
end