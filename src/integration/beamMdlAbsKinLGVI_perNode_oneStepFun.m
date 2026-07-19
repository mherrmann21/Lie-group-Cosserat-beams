function [R_k1, x_k1, eta_k, mu_k, xi_k, solData] = ...
        beamMdlAbsKinLGVI_perNode_oneStepFun( ...
        R_k, x_k, eta_k0, mu_k0, xi_k0, h, simPars, params, isCantilever, ...
        solverConfig)
    % LGVI function to integrate over one timestep (from k -> k1)
    % Contains the integrator from Demoures 2015 without external forces;
    % linear strain-rate damping is included. The solution scheme is the
    % "multi-step" solution proposed in Demoures 2015.
    %
    % Inputs: See below.
    %
    % Outputs:
    %   R_k1    Array of node rotation matrices at next time step
    %   x_k1    Array of node position vectors at next time step
    %   eta_k   Array of discrete node velocities in interval [k, k1]
    %   mu_k    Array of discrete node momenta at current time step
    %   xi_k    Array of discrete deformation gradients at current time
    %           step
    %   solData Struct with solver metadata from the current time step
    %
    % Maximilian Herrmann
    % Chair of Automatic Control
    % TUM School of Engineering and Design
    % Technical University of Munich

    arguments
        % Array of node rotation matrices at current time step k (3, 3, nNodes)
        R_k (3,3,:) double

        % Array of node position vectors at current time step k (3, nNodes)
        x_k (3,:) double

        % Array of discrete velocities at previous time step k0 (6, nNodes)
        eta_k0 (6,:) double

        % Array of discrete momentum mu at previous time step k0 (6, nNodes)
        mu_k0 (6,:) double

        xi_k0 (6,:) double

        % Time step / Sample time
        h (1,1) double

        % Simulation parameters
        simPars (1,1) beamSimPars

        % Parameter object
        params (1,1) beamParams

        % Specifies if the beam is a cantilever, i.e., if first node is
        % fixed
        isCantilever (1,1) logical

        % Struct containing solver configs
        solverConfig (1,1) beamSolverConfig
    end

    % Nr. of nodes and segments
    nNodes = size(simPars.g0, 3);
    nSeg   = nNodes - 1;

    %%% Metadata vectors/matrices
    % Residual error of the implicit solution
    % Number of iterations needed for the implicit solution
    % Exit flag of the implicit solver
    solData.ImplicitError       = nan(nNodes, 1);
    solData.ImplicitIterations  = nan(nNodes, 1);
    solData.ExitFlag            = nan(nNodes, 1);

    % Explicitly define loop variables for code generation

    R_k1 = zeros(3,3,nNodes);
    x_k1 = zeros(3, nNodes);

    mu_k  = zeros(6,nNodes);
    eta_k = zeros(6,nNodes);

    % If the beam is a cantilever: Skip first node and set its velocities
    % / other values to zero / the fixed values
    if isCantilever
        nStart  = 2;

        R_k1(:,:, 1) = simPars.g0(1:3, 1:3,1);
        x_k1  (:, 1) = simPars.g0(1:3, 4,  1);
        mu_k (:, 1) = zeros(6,1);
        eta_k(:, 1) = zeros(6,1);
    else
        nStart = 1;
    end

    % Compute discrete deformations
    % Note: The computed values for xi correspond to discrete updates, not
    % gradients; thus, they are divided by segment length l afterwards
    xi_k = computeDiscreteDeformations( ...
        SE3Matrix(R_k(:,:,:), x_k(:,:))) / (params.L / nSeg);

    % Compute segment length
    l = params.L/ (nNodes - 1);

    % Pre-compute stress and damping force for all nodes (in the "relative"
    % space)
    tauStrDmp = ...
        + params.Cgen * (xi_k - simPars.xiRef) ...
        + diag(params.d)*(xi_k - xi_k0);

    for iN = nStart:nNodes
        
        %%% Strain term: Stiffness and Damping
        if iN == 1
            % First Node
            fStrDmp_i = 2 * cayRTDInvSE3(xi_k(:, iN) * l).' * tauStrDmp(:, iN);
        elseif iN == nNodes
            % Last Node
            fStrDmp_i = - 2 * lAdSE3( caySE3(xi_k(:, iN-1) * l) ).' * ...
                cayRTDInvSE3(xi_k(:, iN-1) * l).' * tauStrDmp(:, iN-1);
        else
            % Interior nodes
            fStrDmp_i = ...
                + cayRTDInvSE3(xi_k(:, iN) * l).' * tauStrDmp(:, iN) ...
                - lAdSE3( caySE3( xi_k(:, iN-1) * l ) ).' ...
                * cayRTDInvSE3(xi_k(:, iN-1) * l).' * tauStrDmp(:, iN-1);
        end

        %%% Update step 1: Compute mu_k from EOM
        mu_k(:,iN) = (...
            + ColAdSE3( caySE3( eta_k0(:,iN)*h ) ) * mu_k0(:,iN) ...
            - h * [ zeros(3,1); R_k(:,:,iN).' * params.m * simPars.g * [0;0;1]] ...
            + h/l * fStrDmp_i ...
            );


        %%% Update step 2: Solve implicitly for eta_k

        tStepFun = @(eta_a_k) beamMdlMu(eta_a_k, h, params) - mu_k(:,iN);
        tStepJac = @(eta_a_k) tStepFunJac_muOnly(eta_a_k, h, params.J, params.m);

        % Initial values for implicit solution: Use eta from last step
        [ ...
            eta_k(:,iN),  ...
            numIterations, ...
            errorFlag, ...
            normFun, ...
            ~, ~ ...
            ] = broydenGood( ...
            tStepFun, tStepJac, eta_k0(:,iN), ...
            solverConfig.errorMargin, solverConfig.maxIterations ...
            );


        %%% Update step 3: Update configuration vector

        % Compute updated SE3 element
        g_i_k1 = SE3Matrix( R_k(:,:,iN), x_k(:,iN) ) * caySE3( eta_k(:,iN)*h);

        % Assign R, x to arrays
        [R_k1(:,:,iN), x_k1(:,iN)] = RxFromSE3Matrix(g_i_k1);


        %%% Housekeeping and Statistics
        solData.ImplicitError(iN)      = normFun;
        solData.ImplicitIterations(iN) = numIterations;
        solData.ExitFlag(iN)           = errorFlag;

    end
end

function mu = beamMdlMu(eta_a_k, h, params)
    % Compute Eta expression for the discrete kinetic energy of a node a
    %
    % Inputs: See below.

    arguments
        % Vector in R6 with the discrete velocity of node a at time step k
        eta_a_k (6,1) double

        % Time step
        h       (1,1) double

        % Parameter object
        params  (1,1) beamParams
    end

    mu = cayRTDInvSE3(eta_a_k*h).' * params.Mgen * eta_a_k;

end
