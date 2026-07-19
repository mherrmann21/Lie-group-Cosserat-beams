function [simData, metaDataSteps] = beamMdlRelKinVarInt_Recursive( ...
        simPars, beamPars, solverConfig, Ba, Bc ) %#codegen
    %% Variational Integrator for the GE beam in relative kinematic description and Reduced Deformation Modes
    %
    % Uses the linear-time recursive variational integrator formulation
    % from [Lee+20] for chain-kinematic multibody systems.
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
    R   = zeros(3,3, nNodes, nSteps+1); % Configuration / R
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

    % Note: Reduced deformations psi don't need to be stored; just keep
    % them in one variable that is updated every iteration
    xi_k  = xi(:,:,1);
    xi_k0 = xi(:,:,1);

    % Compute relative spatial updates for first step
    g_xi_k1 = zeros(4,4,nSeg);
    for iN = 1:nSeg
        g_xi_k1(:,:,iN) = caySE3( xi(:, iN) * discPars.l );
    end


    %% Integration loop

    nStepsDone = nSteps;

    % Initial assignment outside the loop for code generation
    eta_k = eta(:,:,1);

    for k = 1:nSteps

        % Compute external node forces for current time step
        [f_node_k_b, f_node_k_s] = getExtStepNodeForces(simPars, tout(k));

        % Force solver iteration at first time step to avoid convergence
        % problems when the beam is excited from equilibrium
        forceSolverIteration = k == 1;

        [...
            R(:,:,:,k+1), ...
            x(:,:,k+1), ...
            eta_k, ...
            xi_k1, ...
            g_xi_k1, ...
            solData_k ] = beamMdlRelKinVarInt_Recursive_oneStepFun( ...
            R(:,:,:,k), ...
            xi_k, ...
            eta_k, ...
            xi_k0, g_xi_k1, ...
            f_node_k_b, f_node_k_s, ...
            forceSolverIteration, ...
            simPars, discPars, solverConfig, Ba, Bc ...
            );

        eta(:,:,k) = eta_k;
        xi(:,:,k+1) = xi_k1;

        % Assign xi variables to correct time steps for next iteration
        xi_k0 = xi_k;
        xi_k  = xi_k1;

        % Housekeeping and Statistics
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

    % Discrete velocities and momentum at the final step:
    % Not defined since there is no future time step anymore
    % (velocity at k is the velocity in interval k, k+1)
    eta(:,:,end) = nan(6, nNodes);


    %% Assign to output object

    simData = beamSimData;
    simData.R   = R(:,:,:,1:nStepsDone+1);
    simData.x   = x  (:,:,1:nStepsDone+1);
    simData.eta = eta(:,:,1:nStepsDone+1);
    simData.xi  = xi (:,:,1:nStepsDone+1);
    simData.tout = tout(1:nStepsDone+1, 1);

    metaDataSteps = beamSimMetaDataSteps;
    metaDataSteps.ImplicitError      = ImplicitError(1:nStepsDone+1);
    metaDataSteps.ImplicitIterations = ImplicitIterations(1:nStepsDone+1);
    metaDataSteps.ExitFlag           = ExitFlag(1:nStepsDone+1);


end


%% Local Functions for the Beam Model

function [ ...
        R_k1, x_k1, eta_k, xi_k1, g_xi_k1, solData ...
        ] = beamMdlRelKinVarInt_Recursive_oneStepFun( ...
        R_k, xi_k, eta_k0, xi_k0, g_xi_k, ...
        f_node_k_b, f_node_k_s, ...
        forceSolverIteration, ...
        simPars, discPars, solverConfig, ...
        Ba, Bc )
    %% Variational integrator function to integrate over one timestep
    % (from k -> k1), using the recursive algorithms from [Lee+20].
    %
    % Maximilian Herrmann
    % Chair of Automatic Control
    % TUM School of Engineering and Design
    % Technical University of Munich

    arguments
        % Array of node rotation matrices at current time step k (3, 3, nNodes)
        R_k     (3,3,:) double

        % Array of discrete segment deformations (6, nSeg)
        xi_k    (6,:) double

        % Array of discrete velocities at previous time step k0 (6, nNodes)
        eta_k0  (6,:) double

        % Array of allowed discrete segment deformations (6, nSeg)
        % in the last time step
        xi_k0  (:,:) double

        % Array of relative node configurations (in the chain)
        % Index a corresponds to the relative update between node a-1 and a
        g_xi_k  (4,4,:) double

        % External node forces (wrenches in se3*) at current time step (6,nNodes)
        f_node_k_b   (6,:) double   % Wrench in the local / body frame
        f_node_k_s   (6,:) double   % Wrench in the inertial/spatial frame

        % Force a solver iteration by adding a small perturbation to the
        % Initial value of the implicit solver
        forceSolverIteration (1,1) logical

        % beamSimPars object containing the simulation parameters
        simPars (1,1) beamSimPars

        % Parameters of the discrete beam
        discPars (1,1) beamParamsDiscrete

        % Struct containing solver configs
        solverConfig (1,1) beamSolverConfig

        % Selection matrix for allowed discrete deformations (6, nAllwd)
        Ba      (6, :) double

        % Selection matrix for constrained / constant discrete deformations
        % (6, 6-nAllwd)
        Bc      (6, :) double
    end

    % Get variables
    h  = simPars.h;

    % Nr. of nodes, segments and allowed deformation DoFs
    nAllwd = discPars.nAllwd;
    nSeg   = discPars.nSeg;

    % segment length
    l = discPars.l;


    %% Forward Kinematics

    % Array holding all node Jacobians
    % Important: First node is omitted since it is assumed fixed!
    SNode_k = zeros(6, nAllwd, nSeg);

    for ii = 1:nSeg
        % Node Jacobian
        SNode_k(:,:,ii) = l * cayRTDSE3( -xi_k(:, ii) * l ) * Ba;
    end

    % Precompute the right-trivialized derivative for xi_k (if required),
    % so that we don't have to do that in the solver loop
    % Directly multiply with Ba.' to get the reduced matrix
    dTaoInvXi = zeros(nAllwd,6,nSeg);
    if any(discPars.DLinRed(:))
        for iSeg = 1:nSeg
            dTaoInvXi(:,:,iSeg) = Ba .' * cayRTDInvSE3(l*xi_k(:,iSeg));
        end
    end


    %% Node Forces
    % i.e., EOM terms that are multiplied with the Jacobian

    for iN = 2:(nSeg+1)
        f_node_k_b(:, iN) = ...
            ... % Inertia term
            - cayRTDInvSE3(-eta_k0(:,iN)*h).' * discPars.MgenNode(:,:,iN) * eta_k0(:,iN)...
            ...% External body-fixed forces
            - h*f_node_k_b(:, iN) ...
            ... % External spatial forces
            + h*[
            R_k(:,:,iN).' * -f_node_k_s(1:3,iN)
            R_k(:,:,iN).' * -f_node_k_s(4:6,iN)
            ] ...
            ...% Gravity
            + h*simPars.g * [
            discPars.m_a(iN) * cross( discPars.x_a(:,iN), R_k(:,:,iN).' * [0;0;1] )
            discPars.mNode(iN)*R_k(:,:,iN).' * [0;0;1] 
            ];
    end


    %% Segment Forces
    % I.e., terms that are not multiplied with the Jacobian

    % Ext. forces and stress
    f_seg_k = discPars.CgenRed * Ba.' * (xi_k - simPars.xiRef);


    %% Find Root of implicit DEL equation
    % Using the Linear-Time Root Updating algorithm from [Lee+20],
    % Algorithm 3

    % Initial value for the solution
    % Compute as explicit Euler step as done in [Lee+20, Sec.3.3], IG2
    xi_k1 = 2*xi_k - xi_k0;
    psi_k1 = Ba.' * xi_k1;

    % Initialize ExitFlag to 1; will be set to 0 if a solution is found
    solData.ExitFlag = 1;
    solData.ImplicitIterations = solverConfig.maxIterations;

    % Initialize loop variables for code generation
    resNorm = 0;
    eta_k   = zeros(6,nSeg);
    g_xi_k1 = zeros(4,4,nSeg);


    for iIteration = 1:solverConfig.maxIterations

        % Evaluate DEL to get residual
        [resDEL, eta_k,  g_xi_k1, CayRTDMassMatrix_k] = beamMdlRelKinVarInt_DEL_Recursive( ...
            xi_k1, g_xi_k, SNode_k, discPars, h, -f_node_k_b(:,2:end), f_seg_k, dTaoInvXi ...
            );

        % Check residual
        resNorm = norm(resDEL);
        if resNorm > solverConfig.errorMargin || forceSolverIteration
            forceSolverIteration = false;

            % Update estimate of psi_k1 using residual impulse
            % Note: Factor h in the residual impulse is omitted (compared
            % to the algorithm in [Lee+20]) to get fast convergence;
            % maybe this is because the factor h is already in the mass
            % matrix from the ABI algorithm
            psi_k1 = ...
                + psi_k1 ...
                - invMassMatrixTermABI(SNode_k, g_xi_k, resDEL, CayRTDMassMatrix_k, h );

            % Compute full discrete deformations xi from allowed and constrained
            % deformations psiA and psiC
            xi_k1 = Ba * psi_k1 + Bc * Bc.' * simPars.xiRef;
        else
            % Solution found
            % Set iteration count, exit flag and exit loop
            solData.ExitFlag = 0;

            % Note: iteration count is reduced by one; convergence in the
            % first iteration of the loop corresponds to just one function
            % evaluation and no "real" implicit solver iteration
            solData.ImplicitIterations = iIteration-1;
            break;
        end
    end

    % Set implicit error
    solData.ImplicitError = resNorm;

    % Add first node to velocity and momentum arrays (assumed fixed)
    eta_k = [zeros(6,1), eta_k];


    %% Relative Kinematics for the next time step

    % Array of *absolute* node configurations
    % Important: First node is omitted since it is assumed fixed!
    g_k1 = zeros(4,4,nSeg);

    for ii = 1:nSeg
        % Absolute configuration and velocity
        if ii > 1
            g_k1(:,:,ii) = g_k1(:,:,ii-1) * g_xi_k1(:,:, ii);
        else
            % Absolute configuration of the first node includes fixed
            % reference g0
            g_k1(:,:,ii) = simPars.g0(:,:,1) * g_xi_k1(:,:, ii);
        end
    end

    % Explicitly get R and x for current time step
    % Note: In R_k and x_k, the first (fixed) node is *included* for
    % compatibility with the output variables!
    [R_k1, x_k1] = RxFromSE3Matrix(cat(3, simPars.g0(:,:,1), g_k1));
end


function [f_seg_k, eta_k, g_ij_k1, CayRTDMassMatrix_k] = beamMdlRelKinVarInt_DEL_Recursive( ...
        xi_k1, g_ij_k, SNode_k, discPars, h, f_node_k, f_seg_k, dTaoInvXi ...
        )
    %% Recursive Evaluation of the DEL Equations / Discrete Inverse Dynamics
    % Compute the DEL equations of the relative beam model using the
    % Discrete Recursive Euler-Newton Algortithm from [Lee+20] (Algorithm
    % 2)
    %
    % Does not contain the "full" DEL Equs.; only the implicit terms
    % (the explicit terms, that can be computed from known quantities at a
    % time step, must be computed beforehand).
    %
    % Maximilian Herrmann
    % Chair of Automatic Control
    % TUM School of Engineering and Design
    % Technical University of Munich

    arguments (Input)
        % Discrete deformations of all segments at the next time step k+1
        % in the form (6, nSeg)
        xi_k1       (6,:) double

        % Array of SE3 relative node configuration matrices
        % Index a corresponds to the update from the parent node a-1 to
        % the current node a (and index a = tau(xi_a*l) holds)
        g_ij_k      (4,4,:) double

        % Array of joint Jacobians with size (6, nDof, nSeg)
        SNode_k     (6,:,:) double

        % Parameters of the discrete beam
        discPars (1,1) beamParamsDiscrete

        % Time step
        h           (1,1) double

        % Array of size (nAllwd, nNodes) with the (known) external impulses
        % in cartesian space (se*3); i.e., terms that are multiplied with
        % beam Jacobians
        f_node_k     (:,:) double

        % Array of size (nAllwd, nSeg) with the (known) external impulses
        % in joint space(R^nAllwd); i.e., terms that are *not*
        % multiplied with beam Jacobians
        f_seg_k         (:,:) double

        % Array of size (:,6,nSeg) with the precomputed terms
        % Ba.' * tauRTDInv(xi_a_k * l)
        dTaoInvXi      (:,6,:)
    end
    arguments (Output)
        % Array of size (nAllwd, nSeg) with the residual of the DEL
        % equations (segment forces)
        f_seg_k         (:,:) double

        % Array of absolute node velocities with size (6, nSeg)
        % Note: Excludes first node (assumed fixed)
        eta_k       (6,:) double

        % Array of *relative* node configurations (in the chain)
        % Index a corresponds to the relative update between node a-1 and a
        g_ij_k1     (4,4,:) double

        % Array for the product of retraction map RTD and mass matrix
        CayRTDMassMatrix_k    (6,6,:) double
    end

    nSeg = discPars.nSeg;

    %% Pass 1 / Compute Relative Configurations and Velocities

    % Array of *relative* node configurations (in the chain)
    % Index a corresponds to the relative update between node a-1 and a
    g_ij_k1 = zeros(4,4, nSeg);

    % Array of absolute node velocities
    eta_k = zeros(6,nSeg);

    for ii = 1:nSeg
        % Relative node configuration
        % Index a corresponds to the update from the parent node a-1 to
        % the current node a
        g_ij_k1(:,:, ii) = caySE3( xi_k1(:, ii) * discPars.l);

        if ii > 1
            eta_k(:,ii) = cayInvSE3( ...
                g_ij_k(:,:, ii) \ caySE3(h*eta_k(:,ii-1)) * g_ij_k1(:,:, ii) ...
                ) / h;
        else
            % TODO include given "base" velocity if needed
            eta_k(:,ii) = cayInvSE3( ...
                g_ij_k(:,:, ii) \ g_ij_k1(:,:, ii) ) / h;
        end
    end


    %% Compute damping term

    if any(discPars.DLinRed(:) )

        % Array of deformation velocities for all segments
        psi_dot_k = zeros(discPars.nAllwd, nSeg);
        for ii = 1:nSeg
            if ii > 1
                % Segments 2...n
                psi_dot_k(:,ii) = dTaoInvXi(:,:,ii) * (...
                    -eta_k(:,ii-1) + lAdSE3(g_ij_k(:,:,ii)) * eta_k(:,ii) ...
                    )/discPars.l;
            else
                % First segment: Consider that velocity of first node is 0
                psi_dot_k(:,ii) = dTaoInvXi(:,:,ii) * (...
                    lAdSE3(g_ij_k(:,:,ii)) * eta_k(:,ii) ...
                    )/discPars.l;
            end
        end

        % Add linear dissipation (in strain rates) to relative forces
        f_seg_k = f_seg_k + discPars.DLinRed * psi_dot_k;
    end


    %% Pass 2 / Compute forces and residual

    % Array of interaction forces between current body and parent
    F_i = zeros(6,nSeg);

    % Array for the product of retraction map RTD and mass matrix
    CayRTDMassMatrix_k = zeros(6,6,nSeg);

    % Only evaluate for nodes > 1, since Jacobian of the first node is
    % always zero
    for ii = nSeg:(-1):1

        % Explicitly compute the term of the retraction map RTD and mass
        % matrix to be able to re-use it (computation is expensive)
        CayRTDMassMatrix_k(:,:,ii) = cayRTDInvSE3(eta_k(:,ii)*h).' * discPars.MgenNode(:,:,ii+1);

        % Compute interaction forces
        F_i(:,ii) = ( CayRTDMassMatrix_k(:,:,ii) * eta_k(:,ii) - f_node_k(:,ii)  );

        % Add term from child body
        if ii < nSeg
            F_i(:,ii) = ...
                + F_i(:,ii) ...
                + lAdSE3Inv(g_ij_k(:,:, ii + 1)).' * F_i(:,ii + 1);
        end

        % Compute residual
        f_seg_k(:,ii) = SNode_k(:,:,ii).' * F_i(:,ii) + f_seg_k(:,ii)*h;
    end
end
