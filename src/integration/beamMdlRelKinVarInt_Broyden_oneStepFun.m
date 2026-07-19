function [R_k1, x_k1, eta_k, xi_k1, solData, H_k, g_xi_k1, JBeam_k1] = beamMdlRelKinVarInt_Broyden_oneStepFun( ...
        R_k, xi_k, eta_k0, xi_k0, g_xi_k,...
        simPars, H_k, updateInvJacobian, forceSolverIteration, ...
        f_node_k_b, f_node_k_s, f_seg_k, ...
        discPars, solverConfig, ...
        Ba, xiC, JBeam_k ) %#codegen
    %% Variational integrator function to integrate over one timestep
    % (from k -> k1).
    %
    % Inputs: See below.
    %
    % Outputs:
    %   R_k1, x_k1  Node rotations and positions at the next time step
    %   eta_k       Node velocities at the current time step
    %   xi_k1       Segment deformations at the next time step
    %   solData     Solver convergence metadata
    %   H_k         Approximation of the inverse residual Jacobian
    %   g_xi_k1     Relative segment transformations at the next time step
    %   JBeam_k1    Geometric beam Jacobians at the next time step
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

        % Array of discrete segment deformations (6, nSeg) at the previous
        % time step
        xi_k0  (6,:) double

        g_xi_k  (4,4,:) double

        % beamSimPars object containing the simulation parameters
        simPars (1,1) beamSimPars

        % Initial approximation of the inverse Jacobian matrix for the
        % implicit equation system
        H_k    (:,:) double

        % If true, compute the (inverse) Jacobian instead of using the
        % given inverse Jacobian H_k
        updateInvJacobian (1,1) logical

        % Force a solver iteration by ignoring the error threshold
        forceSolverIteration (1,1) logical

        % External node forces (wrenches in se3*) at current time step (6,nNodes)
        f_node_k_b   (6,:) double   % Wrench in the local / body frame
        f_node_k_s   (6,:) double   % Wrench in the inertial/spatial frame

        % Segment forces at current time step (nAllwd,nSeg)
        f_seg_k   (:,:) double

        % Parameters of the discrete beam
        discPars (1,1) beamParamsDiscrete

        % Struct containing solver configs
        solverConfig (1,1) beamSolverConfig

        % Selection matrix for allowed discrete deformations (6, nAllwd)
        Ba      (6, :) double

        % Constant parts of the strain vector (6,nSeg)
        xiC  (6,:) % = Bc .* psiC;

        % Geometric beam Jacobian at current time step
        JBeam_k (:,:,:)
    end

    % Get variables
    h  = simPars.h;

    % Nr. of nodes, segments and allowed deformation DoFs
    nAllwd = discPars.nAllwd;
    nSeg   = discPars.nSeg;
    nNodes = discPars.nNodes;

    % segment length
    l = discPars.l;


    %% Precompute the right-trivialized derivative for xi_k (if required),
    % so that we don't have to do that in the solver loop;
    % Directly multiply with Ba.' to get the reduced matrix
    dTaoInvXi = zeros(nAllwd,6,nSeg);
    if any(discPars.DLinRed(:))
        for iSeg = 1:nSeg
            dTaoInvXi(:,:,iSeg) = Ba .' * cayRTDInvSE3(l*xi_k(:,iSeg));
        end
    end

    %% Node Forces
    % i.e., EOM terms that are multiplied with the Jacobian

    for iN = 2:nNodes
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
    f_seg_k = f_seg_k + discPars.CgenRed * Ba.' * (xi_k - simPars.xiRef);


    %% Solve implicit function

    % Initial value for the solution
    % Compute as explicit Euler step as done in [Lee+20, Sec.3.3], IG2
    xi_k1 = 2*xi_k - xi_k0;
    psi_k1 = zeros(nAllwd,nSeg);% Required for Simulink Codegen
    psi_k1 = Ba.' * xi_k1;

    % Initialize ExitFlag to 1; will be set to 0 if a solution is found
    solData.ExitFlag = 1;
    solData.ImplicitIterations = solverConfig.maxIterations;

    %%% Initial function evaluation

    % Evaluate DEL to get residual
    [resDEL, eta_k, g_xi_k1] = beamMdlRelKinVarInt_Broyden_implicitFun( ...
        xi_k1, JBeam_k, discPars, h, f_node_k_b, f_seg_k(:), g_xi_k, dTaoInvXi ...
        );

    %%% Actual solver loop
    resNorm = norm(resDEL);
    if resNorm > solverConfig.errorMargin || forceSolverIteration

        % Update Implicit Jacobian Matrix if Necessary
        if updateInvJacobian
            % Compute mass matrix and absolute dissipation term
            MBeam = zeros(nAllwd*nSeg);
            for iN = 2:nNodes
                MBeam = MBeam ...
                    + JBeam_k(:,:,iN).' * (...
                    + cayRTDInvSE3( eta_k0(:, iN) * h ).' * discPars.MgenNode(:,:,iN) ...
                    + 2 * h * diag( discPars.dQuad .*abs(eta_k0(:, iN)) ) ...
                    ) * JBeam_k(:,:,iN);
            end
            % Add Jacobian term due to linear strain-rate dissipation and
            % invert matrix
            H_k = inv( MBeam/h + discPars.DLinRedSys );
        end

        for iIteration = 1:solverConfig.maxIterations

            % Variables from last solver iteration
            psi_k1_l0 = psi_k1;
            resDEL_l0 = resDEL;

            %%% Apply state update
            psi_k1 = reshape( psi_k1_l0(:) - H_k * resDEL(:), nAllwd, nSeg);

            %%% Compute new residual

            % Compute full discrete deformations xi from allowed and constrained
            % deformations psiA and psiC
            xi_k1 = Ba * psi_k1 + xiC;

            % Evaluate DEL to get residual
            [resDEL, eta_k, g_xi_k1] = beamMdlRelKinVarInt_Broyden_implicitFun( ...
                xi_k1, JBeam_k, discPars, h, f_node_k_b, f_seg_k(:), g_xi_k, dTaoInvXi ...
                );

            %%% Check residual and update H_k
            resNorm = norm(resDEL);
            if resNorm <= solverConfig.errorMargin
                % Solution found
                % Set iteration count, exit flag and exit loop
                solData.ExitFlag = 0;
                solData.ImplicitIterations = iIteration;
                break;
            else
                % Update approximation of the jacobian
                s_k = psi_k1(:) - psi_k1_l0(:);
                y_k = resDEL(:) - resDEL_l0(:);
                H_k = H_k + ((s_k - H_k*y_k) * (s_k.' * H_k)) / (s_k.' * H_k * y_k);
            end
        end
    else
        % First evaluation already satisfies tolerance
        solData.ExitFlag = 0;
        solData.ImplicitIterations = 0;
    end

    % Set implicit error
    solData.ImplicitError = resNorm;

    % Add first node to velocity array (assumed fixed)
    eta_k = [zeros(6,1), eta_k];

    %% Relative Kinematics for the next time step

    % Beam Jacobians and forward kinematics
    [R_k1, x_k1, ~, JBeam_k1] = beamRelKinReducedComplete(xi_k1, g_xi_k1, simPars.g0(:,:,1), l, Ba);


end

function [f_seg_k, eta_k, g_xi_k1] = beamMdlRelKinVarInt_Broyden_implicitFun( ...
        xi_k1, JBeam, discPars, h, f_node_k, f_seg_k, g_xi_k, dTaoInvXi ...
        )
    % Implicit part of the discrete equations of motion for the relKin
    % variational integrator.
    %
    % Used for both the full and the reduced model! (reduced deformation
    % doFs)
    %
    % Does not contain the "full" discrete eom; only the implicit terms
    % (the explicit terms, that can be computed from known quantities at a
    % time step, must be computed beforehand).
    %
    % Maximilian Herrmann
    % Chair of Automatic Control
    % TUM School of Engineering and Design
    % Technical University of Munich

    arguments (Input)
        % Discrete deformations of all segments at the next time step k+1
        % in the form (6,nSeg)
        xi_k1         (6,:) double

        % Array of Beam Jacobians (6, nAllwd*nSeg, nNodes)
        JBeam         (6, :, :) double

        % Parameters of the discrete beam
        discPars    (1,1)   beamParamsDiscrete

        % Time step
        h             (1,1) double

        % Array of size (nAllwd, nNodes) with the known terms that are
        % multiplied with (transposed) beam Jacobians
        f_node_k         (6,:) double

        % Array of size (nAllwd*nSeg,1) with the known terms that are *not*
        % multiplied with beam Jacobians ("segment/relative/joint" forces)
        f_seg_k           (:,1) double

        % Array of size (4,4,nSeg) with the spatial update matrices
        % cay(xi_a_k * l) at time step k
        g_xi_k        (4,4,:) double

        % Array of size (:,6,nSeg) with the precomputed terms
        % Ba.' * tauRTDInv(xi_a_k * l)
        dTaoInvXi      (:,6,:)
    end
    arguments (Output)
        % Vector of size (nAllwd*nSeg,1) with the values of the discrete
        % equations of motion (residual)
        f_seg_k   (:,1) double

        eta_k   (6,:) double

        % Array of size (4,4,nSeg) with the spatial update matrices
        % cay(xi_a_k * l) at time step k+1
        g_xi_k1
    end

    %% Compute absolute velocities

    nSeg = discPars.nSeg;

    % Array of absolute node velocities and spatial update matrices
    % (starting at node 2 -- node 1 is assumed fixed)
    eta_k = zeros(6,nSeg);
    g_xi_k1 = zeros(4,4,nSeg);

    for ii = 1:nSeg
        g_xi_k1(:,:,ii) = caySE3( xi_k1(:, ii) * discPars.l);
        if ii > 1
            eta_k(:,ii) = cayInvSE3( ...
                g_xi_k(:,:,ii) \ caySE3(h*eta_k(:,ii-1)) * g_xi_k1(:,:,ii) ...
                ) / h;
        else
            % TODO include given "base" velocity if needed
            eta_k(:,ii) = cayInvSE3( ...
                g_xi_k(:,:,ii) \ g_xi_k1(:,:,ii) ) / h;
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
                    -eta_k(:,ii-1) + lAdSE3(g_xi_k(:,:,ii)) * eta_k(:,ii) ...
                    )/discPars.l;
            else
                % First segment: Consider that velocity of first node is 0
                psi_dot_k(:,ii) = dTaoInvXi(:,:,ii) * (...
                    lAdSE3(g_xi_k(:,:,ii)) * eta_k(:,ii) ...
                    )/discPars.l;
            end
        end

        % Add linear dissipation (in strain rates) to relative forces
        f_seg_k = f_seg_k + discPars.DLinRedSys * psi_dot_k(:);
    end


    %% Evaluate DEL / Compute Residual


    f_seg_k = f_seg_k*h;
    % Only evaluate for nodes > 1, since Jacobian of the first node is
    % always zero
    for iN = 2:discPars.nNodes

        f_seg_k = f_seg_k ...
            + JBeam(:, :, iN).' *( ...
            + cayRTDInvSE3(eta_k(:,iN-1)*h).' * discPars.MgenNode(:,:,iN) * eta_k(:,iN-1) ...
            + f_node_k(:, iN) ...
            ... % Quadratic dissipation in absolute velocities
            + h* discPars.dQuad .* eta_k(:,iN-1).^2 .* sign(eta_k(:,iN-1)) ...
            );
    end
end
