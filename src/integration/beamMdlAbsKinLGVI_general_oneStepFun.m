function [R_k1, x_k1, eta_k, xi_k1, solData, H_k] = ...
        beamMdlAbsKinLGVI_general_oneStepFun( ...
        R_k, x_k, eta_k0, xi_k, simPars, ...
        H_k, updateInvJacobian, forceSolverIteration, params, ...
        f_node_k_b, f_node_k_s, ...
        solverConfig)
    % LGVI function to integrate over one timestep (from k -> k1)
    % Contains the integrator from Demourez 2015 with external forces and
    % linear strain-rate dissipation;
    % the DEL-equs. are solved directly with Broyden's good method.
    %
    % Inputs: See below.
    %
    % Outputs:
    %   R_k1    Array of node rotation matrices at next time step
    %   x_k1    Array of node position vectors at next time step
    %   eta_k   Array of discrete node velocities in interval [k, k1]
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
        R_k     (3,3,:) double

        % Array of node position vectors at current time step k (3, nNodes)
        x_k     (3,:) double

        % Array of discrete velocities at previous time step k0 (6, nNodes)
        eta_k0  (6,:) double

        % Array of discrete deformations at current time step k (6,nSeg)
        xi_k    (6,:) double

        % beamSimPars object containing the simulation parameters
        simPars (1,1) beamSimPars

        % Initial approximation of the inverse Jacobian matrix for the
        % implicit equation system
        H_k    (:,:) double

        % If true, compute the (inverse) Jacobian instead of using the
        % given inverse Jacobian H_k
        updateInvJacobian (1,1) logical

        % Force a solver iteration by adding a small perturbation to the
        % Initial value of the implicit solver
        forceSolverIteration (1,1) logical

        % Parameter object
        params  (1,1) beamParams

        % Specifies if the beam is a cantilever, i.e., if first node is
        % fixed
        %isCantilever (1,1) logical

        % External node forces (wrenches in se3*) at current time step (6,nNodes)
        f_node_k_b   (6,:) double   % Wrench in the local / body frame
        f_node_k_s   (6,:) double   % Wrench in the inertial/spatial frame

        % Struct containing solver configs
        solverConfig (1,1) beamSolverConfig
    end

    % Get variables
    h  = simPars.h;

    % Nr. of nodes and segments
    nNodes = size(simPars.g0, 3);
    nSeg   = nNodes - 1;

    % If the beam is a cantilever: Skip first node and set its velocities
    % / other values to zero / the fixed values
    % Currently, we only consider cantilevers.
    isCantilever = 1;
    nStart = 2;
    % if isCantilever
    %     nStart = 2;
    % else
    %     nStart = 1;
    % end

    % Compute segment length
    l = params.L / nSeg;


    %% Compute known terms

    % Node forces ("absolute" forces)

    g_k = SE3Matrix( R_k, x_k);

    fac = [2; ones(nNodes-2, 1); 2]*h/l;
    for iN = nStart:nNodes
        f_node_k_b(:, iN) = ...
            ...
            - fac(iN) * f_node_k_b(:, iN) ...
            - cayRTDInvSE3(-eta_k0(:,iN)*h).' * params.Mgen * eta_k0(:,iN)...
            + [
            R_k(:,:,iN).' *  - fac(iN) * f_node_k_s(1:3,iN);
            R_k(:,:,iN).' * (- fac(iN) * f_node_k_s(4:6,iN) ...
            + h * params.m * simPars.g * [0;0;1] )
            ];
    end

    % Combined strain terms in the DEL equations
    Q_k = zeros(6,nNodes);

    % Precompute matrix VUDUV for dissipation terms
    % Last index:
    %    1 = block on the diagonal
    %    2 = block to the right of diag. = transposed block below diag.
    D   = diag(params.d);
    UV  = zeros(6,6,2);
    UV0 = zeros(6,6,2);
    VUDUV = zeros(6,6,nNodes,2);

    for iSeg = 1:nSeg
        cay = cayRTDInvSE3( l*xi_k(:,iSeg) );
        UV(:,:,1) = -cay;
        UV(:,:,2) = cay * lAdSE3( caySE3( l*xi_k(:, iSeg) ));
        if iSeg == 1
            % First node
            Q_k(:,iSeg)  = UV(:,:,1).' * params.Cgen * (xi_k(:,iSeg) - simPars.xiRef(:,iSeg));
            VUDUV(:,:,1) = UV(:,:,1).' * D * UV(:,:,iSeg,1);
        else
            % Inner Nodes
            Q_k(:,iSeg) = ...
                + UV0(:,:,2).' * params.Cgen * (xi_k(:,iSeg-1) - simPars.xiRef(:,iSeg-1))...
                +  UV(:,:,1).' * params.Cgen * (xi_k(:,iSeg)   - simPars.xiRef(:,iSeg));
            VUDUV(:,:,iSeg,1) = ...
                + UV0(:,:,2).' * D * UV0(:,:,2)...
                +  UV(:,:,1).' * D * UV(:,:,1);
        end
        VUDUV(:,:,iSeg,2) = UV(:,:,1).' * D * UV(:,:,2);
        UV0 = UV;
    end

    % Last node
    Q_k(:,nNodes) = UV(:,:,2).' * params.Cgen * (xi_k(:,nSeg) - simPars.xiRef(:,nSeg));
    VUDUV(:,:,nNodes,1) = UV(:,:,2).' * D * UV(:,:,2);


    %% Solve implicit function

    % Initial value for the solution
    % Use value from last time step
    % TODO: Use an explicit Euler step at the velocity level as done in
    % [Lee+20, Sec.3.3], IG2?
    eta_k = eta_k0;

    % Initialize ExitFlag to 1; will be set to 0 if a solution is found
    solData.ExitFlag = 1;
    solData.ImplicitIterations = solverConfig.maxIterations;

    %%% Initial function evaluation
    resDEL = beamMdlAbsKinDEL( ...
        eta_k, f_node_k_b, Q_k, g_k, VUDUV, h, params, nStart ...
        );

    %%% Actual solver loop

    resNorm = norm(resDEL);
    if resNorm > solverConfig.errorMargin || forceSolverIteration

        % Update Implicit Jacobian Matrix if Necessary
        if updateInvJacobian
            H_k = invDELJac(eta_k, VUDUV, params, h, nStart, solverConfig.UseExactJacobian);
        end

        for iIteration = 1:solverConfig.maxIterations

            % Variables from last solver iteration
            eta_k_l0  = eta_k;
            resDEL_l0 = resDEL;

            % Apply state update
            eta_k = eta_k_l0 - reshape( H_k * resDEL , 6, nNodes);

            if isCantilever
                eta_k(:,1) = zeros(6,1);
            end

            %%% Compute new residual
            resDEL = beamMdlAbsKinDEL( ...
                eta_k, f_node_k_b, Q_k, g_k, VUDUV, h, params, nStart ...
                );

            %%% Check residual
            resNorm = norm(resDEL);
            if resNorm > solverConfig.errorMargin
                % Update approximation of the jacobian
                s_k = eta_k(:)  - eta_k_l0(:);
                y_k = resDEL - resDEL_l0;
                H_k = H_k + ((s_k - H_k*y_k) * (s_k.' * H_k)) / (s_k.' * H_k * y_k);
            else
                % Solution found
                % Set iteration count, exit flag and exit loop
                solData.ExitFlag = 0;
                solData.ImplicitIterations = iIteration;
                break;
            end
        end
    else
        % First evaluation already satisfies tolerance
        solData.ExitFlag = 0;
        solData.ImplicitIterations = 0;
    end

    % Set implicit error
    solData.ImplicitError = resNorm;

    % Compute updated configuration and deformation gradients from discrete
    % velocity
    xi_k1 = zeros(6,nSeg);
    g_k1  = zeros(4,4,nNodes);
    if isCantilever
        g_k1(:,:, 1) = simPars.g0(:,:,1);
    end
    for iN = nStart:nNodes
        % Compute updated SE3 element
        g_k1(:,:,iN) = g_k(:,:,iN) * caySE3( eta_k(:,iN)*h);
    end
    for iN = 1:nSeg
        % Compute discrete deformation gradient
        xi_k1(:,iN) = cayInvSE3(g_k1(:,:,iN) \ g_k1(:,:,iN+1))/l;
    end

    % Assign R, x to arrays
    [R_k1, x_k1] = RxFromSE3Matrix(g_k1);

end

function residual = beamMdlAbsKinDEL(eta_k, F_k0, Q_k, g_k, VUDUV, h, params, nStart)
    %% Evaluate complete DEL equations for the absKin beam Model

    nNodes = size(g_k, 3);
    nSeg = nNodes - 1;
    l = params.L / nSeg;

    residual = zeros(6*nNodes,1);
    fac = [2; ones(nNodes-2, 1); 2]*h/l;
    for iN = nStart:nNodes

        % Compute dissipation terms
        if any(params.d)
            if iN == 1
                d_a = ...
                    + VUDUV(:,:,iN,1) * eta_k(:,iN) ...
                    + VUDUV(:,:,iN,2) * eta_k(:,iN+1);
            elseif iN == nNodes
                d_a = ...
                    + VUDUV(:,:,iN-1,2).' * eta_k(:,iN-1) ...
                    + VUDUV(:,:,iN,1) * eta_k(:,iN);
            else
                d_a = ...
                    + VUDUV(:,:,iN-1,2).' * eta_k(:,iN-1) ...
                    + VUDUV(:,:,iN,1) * eta_k(:,iN) ...
                    + VUDUV(:,:,iN,2) * eta_k(:,iN+1);
            end
        else
            d_a = zeros(6,1);
        end


        % Compute residual
        residual(iN*6-5 : iN*6) = ...
            + cayRTDInvSE3(eta_k(:,iN)*h).' * params.Mgen * eta_k(:,iN) ...
            + F_k0(:,iN) + (d_a/l + Q_k(:,iN))*fac(iN);
    end
end

function H_k = invDELJac(eta_k, VUDUV, params, h, nStart, useExactJacobian)
    % Compute the inverse of the DEL Jacobian
    nNodes = size(eta_k, 2);
    nSeg = nNodes - 1;

    l = params.L / nSeg;
    fac = [2; ones(nNodes-2, 1); 2]*h/l^2;

    % Check if damping is present to decide whether we have to use the
    % full, complicated Jacobian or only have to compute the inertia terms
    if any(params.d)  && useExactJacobian
        J_k = zeros(6,6,nNodes,nNodes);
        for iSeg = 1:nSeg
            J_k(:,:,iSeg,iSeg)   = ...
                + VUDUV(:,:,iSeg,1)* fac(iSeg) ...
                + tStepFunJac_muOnly(eta_k(:,iSeg), h, params.J, params.m);
            J_k(:,:,iSeg,iSeg+1) = VUDUV(:,:,iSeg,2)   * fac(iSeg);
            J_k(:,:,iSeg+1,iSeg) = VUDUV(:,:,iSeg,2).' * fac(iSeg+1);
        end
        J_k(:,:,nNodes,nNodes) = ...
            + VUDUV(:,:,nNodes,1)*fac(nNodes) ...
            + tStepFunJac_muOnly(eta_k(:,nNodes), h, params.J, params.m);

        H_k = inv( reshape(permute(J_k,[1,3,2,4]),[6*nNodes,6*nNodes]) );

    else
        H_k = zeros(6,6,nNodes,nNodes);
        for iN = nStart:nNodes
            H_k(:,:,iN,iN) = inv(...
                tStepFunJac_muOnly(eta_k(:,iN), h, params.J, params.m) ...
                );
        end
        H_k = reshape(permute(H_k,[1,3,2,4]),[6*nNodes,6*nNodes]);
    end
end
