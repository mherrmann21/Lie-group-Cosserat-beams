function [simData, metaDataSteps] = beamMdlAbsKinLGVI_perNode( ...
        simPars, params, solverConfig, Ba, Bc)
    %% LGVI for the GE beam in absolute, global description
    %
    % Inputs: See arguments block below.
    % Outputs: simRes object.
    %
    % Maximilian Herrmann
    % Chair of Automatic Control
    % TUM School of Engineering and Design
    % Technical University of Munich
    %
    % ToDo:
    % * Initial Velocities (not present yet).
    % * External forces

    arguments
        % simPars object (struct) with simulation parameters.
        % See class definition for details.
        simPars         (1,1) beamSimPars

        % Beam parameter struct
        params          (1,1) beamParams

        % Specifies if the beam is a cantilever, i.e., if first node is
        % fixed
        %isCantilever    (1,1) logical

        % Struct containing solver configs
        solverConfig    (1,1) beamSolverConfig

        % Selection matrix for allowed discrete deformations (6, nAllwd)
        Ba   (6, :) double

        % Selection matrix for constrained / constant discrete deformations
        % (6, 6-nAllwd)
        Bc   (6, :) double
    end


    %% Assign input data

    % Time step / Sample time
    h = simPars.h;

    % Nr. of nodes and segments
    nNodes = size(simPars.g0, 3);
    nSeg   = nNodes - 1;

    % Time span [t_0, t_end] (also vector of sample times possible)
    tspan = [0, simPars.tEnd];

    % Nr. of integration steps
    nSteps = round( ( tspan(end) - tspan(1) ) / h );

    % Currently, we only consider cantilevers.
    isCantilever = 1;

    
    %% Initialize output arrays
    % Note: time dimension (outer loop) should be last index,
    % node dimension (inner loop) should be first;
    % Also put data dimensions first, which make squeeze/reshape
    % unnecessary
    R   = zeros(3,3, nNodes, nSteps+1); % Configuration / R
    x   = zeros(3, nNodes,   nSteps+1);   % Configuration / x
    eta = zeros(6, nNodes,   nSteps+1);   % Discrete velocity
    mu  = zeros(6, nNodes,   nSteps+1);   % Momentum
    xi  = zeros(6, nSeg,     nSteps+1);   % Discrete deformation
    
    % Metadata vectors/matrices
    ImplicitError       = nan(nNodes, nSteps+1);
    ImplicitIterations  = nan(nNodes, nSteps+1);
    ExitFlag            = nan(nNodes, nSteps+1);


    %% Assign initial conditions
    R(:,:,:,1) = simPars.g0(1:3, 1:3, :);
    x(:,:,1)   = simPars.g0(1:3, 4, :);
    R(:,:,:,2) = simPars.g0(1:3, 1:3, :);
    x(:,:,2)   = simPars.g0(1:3, 4, :);

    % Compute discrete deformations for initial time step
    % Note: The computed values for xi correspond to discrete updates, not
    % gradients; thus, they are divided by segment length l afterwards
    xi(:,:,1) = computeDiscreteDeformations(simPars.g0) / (params.L / nSeg);
    xi(:,:,2) = xi(:,:,1);

    % Note: Discrete Momentum is (currently) always zero since we assume
    % zero initial velocity


    %% Integration loop

    nStepsDone = nSteps;

    for k = 2:nSteps

        % Function Inputs:
        % (R_k, x_k, eta_k0, mu_k0, xiRef, R0, x0, h, params, isCantilever)
        % Function outputs:
        % [R_k1, x_k1, eta_k, mu_k, xi_k, solData]
        [...
            R(:,:,:,k+1), ...
            x(:,:,k+1), ...
            eta(:,:,k), ...
            mu(:,:,k), ...
            xi(:,:,k), ...
            solData_k...
            ] = beamMdlAbsKinLGVI_perNode_oneStepFun(...
            R(:,:,:,k), x(:,:,k), ...
            eta(:,:,k-1), mu(:,:,k-1), xi(:,:,k-1), ...
            h, simPars, params, isCantilever, solverConfig);

        %%% Housekeeping and Statistics
        ImplicitError(:,k)      = solData_k.ImplicitError;
        ImplicitIterations(:,k) = solData_k.ImplicitIterations;
        ExitFlag(:,k)           = solData_k.ExitFlag;

        % Check if solver was sucessful; cancel simulation if not
        if ~all(isnan(solData_k.ExitFlag)) && any(solData_k.ExitFlag)
            nStepsDone = k;
            break;
        end
    end


    %% Compute variables at final time step (nStep + 1)

    % Compute discrete deformations
    % Note: The computed values for xi correspond to discrete updates, not
    % gradients; thus, they are divided by segment length l afterwards
    xi(:,:,nStepsDone+1) = computeDiscreteDeformations( ...
        SE3Matrix(R(:,:,:,nStepsDone+1), x(:,:,nStepsDone+1)) ...
        ) / (params.L / nSeg);

    % Discrete velocities: Not defined since there is no future time step
    % anymore (velocity at k is the velocity in interval k, k+1)
    eta(:,:,nStepsDone+1) = nan(6, nNodes);

    % Same for discrete momentum
    mu(:,:,nStepsDone+1) = nan(6, nNodes);


    %% Assign to output object

    simData = beamSimData;
    simData.R   = R(:,:,:,1:nStepsDone+1);
    simData.x   = x  (:,:,1:nStepsDone+1);
    simData.eta = eta(:,:,1:nStepsDone+1);
    simData.mu  = mu (:,:,1:nStepsDone+1);
    simData.xi  = xi (:,:,1:nStepsDone+1);
    simData.tout = (tspan(1):h:(tspan(1) + h*nStepsDone) )';

    metaDataSteps = beamSimMetaDataSteps;
    metaDataSteps.ImplicitError      = ImplicitError(:,1:nStepsDone+1);
    metaDataSteps.ImplicitIterations = ImplicitIterations(:,1:nStepsDone+1);
    metaDataSteps.ExitFlag           = ExitFlag(:,1:nStepsDone+1);

end
