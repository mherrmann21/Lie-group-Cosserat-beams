function [simData, metaDataSteps] = beamMdlAbsKinLGVI_general( ...
        simPars, params, solverConfig, Ba, Bc) %#codegen
    %% LGVI for the GE beam in absolute, global description
    %
    % Inputs: See arguments block below.
    %
    % Maximilian Herrmann
    % Chair of Automatic Control
    % TUM School of Engineering and Design
    % Technical University of Munich
    %
    % ToDo:
    % * Initial Velocities (not present yet).

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
        % (Not used)
        Ba   (6, :) double

        % Selection matrix for constrained / constant discrete deformations
        % (6, 6-nAllwd)
        % (Not used)
        Bc   (6, :) double
    end



    %% Assign input data

    % Time step / Sample time
    h = simPars.h;

    % Nr. of nodes and segments
    nNodes = size(simPars.g0, 3);
    nSeg   = nNodes - 1;

    % Nr. of integration steps
    nSteps = round( simPars.tEnd / h );

    % Currently, we only consider cantilevers.
    %isCantilever = 1;


    %% Initialize output arrays
    % Note: time dimension (outer loop) should be last index,
    % node dimension (inner loop) should be first;
    % Also put data dimensions first, which make squeeze/reshape
    % unnecessary
    R   = zeros(3,3, nNodes, nSteps+1); % Configuration / R
    x   = zeros(3,   nNodes, nSteps+1); % Configuration / x
    eta = zeros(6,   nNodes, nSteps+1); % Discrete velocity
    xi  = zeros(6,   nSeg,   nSteps+1); % Discrete deformation

    % Metadata vectors/matrices
    ImplicitError       = nan(1,nSteps+1);
    ImplicitIterations  = nan(1,nSteps+1);
    ExitFlag            = nan(1,nSteps+1);

    % Precompute time vector
    tout = (0:h:h*nSteps )';


    %% Assign initial conditions
    R(:,:,:,1) = simPars.g0(1:3, 1:3, :);
    x(:,:,1)   = simPars.g0(1:3, 4, :);
    R(:,:,:,2) = simPars.g0(1:3, 1:3, :);
    x(:,:,2)   = simPars.g0(1:3, 4, :);

    % Compute discrete deformations for initial time step
    xi(:,:,1) = computeDiscreteDeformations(simPars.g0) / (params.L / nSeg);
    xi(:,:,2) = xi(:,:,1);

    % Initialize variable for the inverse Jacobian matrix of the implicit
    % equations
    % (The value assigned here is not actually used; the inverse Jacobian
    % is computed inside the one-step function for the first time step)
    H_k = zeros(nNodes*6);

    
    %% Integration loop

    nStepsDone = nSteps;

    % Initial assignment outside the loop for code generation
    eta_k = eta(:,:,1);

    for k = 1:nSteps
        if k == 1
            % Initial computation of the Jacobian
            updateInvJacobian = true;

            % Force solver iteration at first time step to avoid convergence
            % problems when the beam is excited from equilibrium
            forceSolverIteration = true;
        else
            forceSolverIteration = false;

            % Check nr. of iterations in the previous time step and recompute
            % Jacobian if necessary
            updateInvJacobian = ImplicitIterations(k-1) > solverConfig.JacobianIterationThreshold;
        end

        % Compute external node forces for current time step
        [f_node_k_b, f_node_k_s] = getExtStepNodeForces(simPars, tout(k));

        % Integrate time step
        [...
            R(:,:,:,k+1), ...
            x(:,:,k+1), ...
            eta_k, ...
            xi(:,:,k+1), ...
            solData_k, ...
            H_k...
            ] = beamMdlAbsKinLGVI_general_oneStepFun(...
            R(:,:,:,k), x(:,:,k), ...
            eta_k, xi(:,:,k), ...
            simPars, H_k, updateInvJacobian, forceSolverIteration, params, ...
            f_node_k_b, f_node_k_s, ...
            solverConfig);

        eta(:,:,k) = eta_k;

        %%% Housekeeping and Statistics
        ImplicitError(k)      = solData_k.ImplicitError;
        ImplicitIterations(k) = solData_k.ImplicitIterations;
        ExitFlag(k)           = solData_k.ExitFlag;

        % Check if solver was successful; cancel simulation if residual is
        % above residual limit
        if ( ~isnan(solData_k.ExitFlag) && solData_k.ExitFlag && ...
                solData_k.ImplicitError > solverConfig.errorMarginLimit ) ...
            || isnan(solData_k.ImplicitError)
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


    %% Assign to output object

    simData = beamSimData;
    simData.R    = R(:,:,:,1:nStepsDone+1);
    simData.x    = x  (:,:,1:nStepsDone+1);
    simData.eta  = eta(:,:,1:nStepsDone+1);
    simData.xi   = xi (:,:,1:nStepsDone+1);
    simData.tout = tout(1:nStepsDone+1, 1);

    metaDataSteps = beamSimMetaDataSteps;
    metaDataSteps.ImplicitError      = ImplicitError(1:nStepsDone+1);
    metaDataSteps.ImplicitIterations = ImplicitIterations(1:nStepsDone+1);
    metaDataSteps.ExitFlag           = ExitFlag(1:nStepsDone+1);
end
