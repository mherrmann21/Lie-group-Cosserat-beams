function [simData, metaDataSteps] = beamMdlRelKinVarInt_Broyden( ...
        simPars, beamPars, solverConfig, Ba, Bc ) %#codegen
    %% Variational Integrator for the GE beam in relative kinematic formulation
    %
    % Maximilian Herrmann
    % Chair of Automatic Control
    % TUM School of Engineering and Design
    % Technical University of Munich

    arguments
        % simPars object (struct) with simulation parameters.
        % See class definition for details.
        simPars         (1,1) beamSimPars

        % Beam parameter struct
        beamPars (1,1) beamParams

        % Struct containing solver configs
        solverConfig (1,1) beamSolverConfig

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
    nAllwd = size(Ba, 2);

    % Nr. of integration steps
    nSteps = round( simPars.tEnd / h );

    % time vector (has length nSteps + 1)
    tout = (0:h:h*nSteps)';


    %% Compute discrete node variables
    discPars = beamParamsDiscrete(simPars, beamPars, Ba);


    %% Initialize output arrays
    % Note: time dimension (outer loop) should be last index,
    % node dimension (inner loop) should be first;
    % Also put data dimensions first, which makes squeeze/reshape
    % unnecessary
    R   = zeros(3,3, nNodes, nSteps+1);   % Configuration / R
    x   = zeros(3, nNodes,   nSteps+1);   % Configuration / x
    eta = zeros(6, nNodes,   nSteps+1);   % Discrete velocity
    xi  = zeros(6, nSeg,     nSteps+1);   % Discrete deformation

    % Metadata vectors/matrices
    ImplicitError       = nan(1,nSteps+1);
    ImplicitIterations  = nan(1,nSteps+1);
    ExitFlag            = nan(1,nSteps+1);


    %% Assign initial conditions
    R(:,:,:,1) = simPars.g0(1:3, 1:3, :);
    x(:,:,1)   = simPars.g0(1:3, 4, :);

    % Constrained strains
    xiC  = Bc * Bc.' * simPars.xiRef;

    % Compute initial strains to be consistent with beam constraints
    xi(:,:,1)  = Ba*Ba.' * computeDiscreteDeformations(simPars.g0) / discPars.l + xiC;

    % Compute beam Jacobian and relative spatial updates for first step
    [JBeam_k1, g_xi_k1] = computeBeamJacobiansReduced(xi(:,:,1), discPars.l, Ba);

    % Initialize variable for the inverse Jacobian matrix of the implicit
    % equations
    % (The value assigned here is not actually used; the inverse Jacobian
    % is computed inside the one-step function for the first time step)
    H_k = zeros(nSeg*nAllwd);


    %% Prepare external forces
    f_seg_k    = zeros(nAllwd,nSeg);

    
    %% Integration loop

    nStepsDone = nSteps;

    % Initial assignments outside the loop for code generation
    eta_k = zeros(6,nNodes);
    xi_k  = xi(:,:,1);
    xi_k0 = xi_k; 
    xi_k1 = xi_k; 

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

            xi_k0 = xi_k;
            xi_k  = xi_k1;
        end

        % Compute external node forces for current time step
        [f_node_k_b, f_node_k_s] = getExtStepNodeForces(simPars, tout(k));

        [...
            R(:,:,:,k+1), ...
            x(:,:,k+1), ...
            eta_k, ...
            xi_k1, ...
            solData_k, ...
            H_k,...
            g_xi_k1, JBeam_k1] = beamMdlRelKinVarInt_Broyden_oneStepFun( ...
            R(:,:,:,k), xi_k,...
            eta_k, xi_k0, g_xi_k1,...
            simPars, H_k, updateInvJacobian, forceSolverIteration, ...
            f_node_k_b, f_node_k_s, f_seg_k, ...
            discPars, solverConfig, ...
            Ba, xiC, JBeam_k1 ...
            );

        eta(:,:,k) = eta_k;
        xi(:,:,k+1) = xi_k1;

        % Housekeeping and Statistics
        ImplicitError(k)      = solData_k.ImplicitError;
        ImplicitIterations(k) = solData_k.ImplicitIterations;
        ExitFlag(k)           = solData_k.ExitFlag;

        % Check if solver was successful; cancel simulation if residual is
        % above residual limit
        if ( solData_k.ExitFlag && solData_k.ImplicitError > solverConfig.errorMarginLimit ) ...
                || isnan(solData_k.ImplicitError)
            nStepsDone = k;
            break;
        end
    end

    % Discrete velocities and momentum at the final step:
    % Not defined since there is no future time step anymore
    % (velocity at k is the velocity in interval k, k+1)
    eta(:,:,end) = nan(6, nNodes);


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
