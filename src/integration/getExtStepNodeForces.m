function [f_node_k_b, f_node_k_s] = getExtStepNodeForces(simPars, t_k)
    %% Get external node forces for a beam for current time step
    %
    % Maximilian Herrmann
    % Chair of Automatic Control
    % TUM School of Engineering and Design
    % Technical University of Munich

    arguments (Input)
        simPars (1,1) beamSimPars

        % Current simulation time
        t_k     (1,1) double
    end
    arguments (Output)
        % Array of external node forces (wrenches) in the body-fixed frame
        % with dimensions (6,nNodes)
        f_node_k_b (6,:) double

        % Array of external node forces (wrenches) in the spatial frame
        % with dimensions (6,nNodes)
        f_node_k_s (6,:) double
    end


    %% Get variables
    nNodes = size(simPars.g0, 3);

    % Body-fixed wrenches
    f_node_k_b = zeros(6, nNodes);

    %% Spatial wrenches
    f_node_k_s = zeros(6, nNodes);

    if ~isempty(simPars.f_node_s)
        switch simPars.force_scaling_mode
            case 0
                % Constant
                f_node_k_s = simPars.f_node_s;
            case 1
                % Smooth impulse
                if t_k <= simPars.force_tEnd
                    f_node_k_s = simPars.f_node_s * ...
                        (1-cos(2*pi*t_k / simPars.force_tEnd ))/2;
                end
            case 2
                % Smooth increase
                if t_k <= simPars.force_tEnd
                    f_node_k_s = simPars.f_node_s * ...
                        (1-cos(pi*t_k / simPars.force_tEnd ))/2;
                else
                    f_node_k_s = simPars.f_node_s;
                end
            case 3
                % Smooth decrease
                if t_k <= simPars.force_tEnd
                    f_node_k_s = simPars.f_node_s * ...
                       (cos(pi*t_k / simPars.force_tEnd)+1)/2;
                end
            case 4
                % Arbitrary interpolation
                if ~isempty(simPars.f_node_tVec) && ~isempty(simPars.f_node_sVec)
                    % If time-dependent interpolation data is given, compute
                    % scaling factor for current time instance and scale
                    % accordingly
                    f_node_k_s = simPars.f_node_s * interp1( ...
                        simPars.f_node_tVec, simPars.f_node_sVec, t_k, ...
                        "linear", 0);
                end
        end
    end

