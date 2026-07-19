function [gEqu, simData, solInfo] = computeStaticEquilibriumRelKin(simPars, beamPars, Ba, Bc, nLoadSteps, useExp)
    %% Compute static equilbrium of a beam
    % Based on the relKin (Reduced) formulation
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
        beamPars          (1,1) beamParams

        % Selection matrix for allowed discrete deformations (6, nAllwd)
        Ba   (6, :) double

        % Selection matrix for constrained / constant discrete deformations
        % (6, 6-nAllwd)
        Bc   (6, :) double

        % Nr. of loading steps
        nLoadSteps (1,1) double

        % Used retraction map: 0 = Cayley, 1 = Exponential Map
        useExp (1,1) logical
    end

    disp('Computing Static Equilibrium Configuration...');


    %% Assign input data
    nAllwd = size(Ba, 2);
    nSeg = size(simPars.xiRef, 2);
    psiC = Bc.' * simPars.xiRef;


    %% Compute discrete node variables
    discPars = beamParamsDiscrete(simPars, beamPars, Ba);


    %% Initialize output arrays
    gEqu = zeros(4,4,nSeg+1);
    REqu = zeros(3,3,nSeg+1);
    xEqu = zeros(3,nSeg+1);
    % xiEqu = zeros(6,nSeg);

    % Metadata vectors/matrices
    solInfo.iterations = nan(nLoadSteps, 1);
    solInfo.exitFlag   = nan(nLoadSteps, 1);
    solInfo.residual   = nan(nLoadSteps, 1);


    %% Assign initial conditions
    xi_k = simPars.xi0;
    psi_k = Ba.' * simPars.xi0;


    %% Solve statics

    opts = optimoptions('fsolve', ...
        'Display','off', ...
        'Algorithm','levenberg-marquardt', ...
        'FunctionTolerance', 1e-14);

    loadFactors = linspace(0,1,nLoadSteps);

    % Solve for equilibrium: Increase load slowly via load steps
    for iLoadStep = 1:nLoadSteps
        if nLoadSteps <9
            fprintf( ...
                '     Starting loading step %.1f/%.1f (factor %.2f)...\n', ...
                iLoadStep, nLoadSteps, loadFactors(iLoadStep) ...
                );
        end

        % Function to solve
        simParsStep = simPars;
        simParsStep.h = 1;
        simParsStep.f_node_b = loadFactors(iLoadStep) * simPars.f_node_b;
        simParsStep.f_node_s = loadFactors(iLoadStep) * simPars.f_node_s;

        simParsStep.g = loadFactors(iLoadStep) * simPars.g;

        equFun = @(psi_equ_Vec) implicitEquilibriumFunction( ...
            reshape(psi_equ_Vec, nAllwd, nSeg), simParsStep, discPars, Ba, Bc, useExp);

        % Solve function
        tic
        [psi_kVec, equResidual, equExitFlag, equMeta] = fsolve(equFun, psi_k(:), opts);
        equTime = toc;

        psi_k = reshape(psi_kVec, nAllwd, nSeg);

        % Compute full discrete deformations xi from allowed and constrained
        % deformations psiA and psiC
        xi_k = Ba * psi_k + Bc * psiC;

        % Compute forward kinematics
        if useExp
            [REqu, xEqu, gEqu, ~] = beamRelFwdKinExp(xi_k *discPars.l, simPars.g0(:,:,1));
        else
            [REqu, xEqu, gEqu, ~] = beamRelFwdKin(xi_k *discPars.l, simPars.g0(:,:,1));
        end

        solInfo.iterations(iLoadStep) = equMeta.iterations;
        solInfo.exitFlag(iLoadStep)   = equExitFlag;
        solInfo.residual(iLoadStep)   = mean(equResidual(:));

        if equExitFlag > 0
            if nLoadSteps <9
                fprintf('     Equilibrium configuration solved in %f s. Iterations: %.1f,  Mean Residual: %e.\n', equTime, equMeta.iterations, mean(equResidual(:)));
            end
        else
            fprintf('     Failed to solve equilibrium configuration.\n');
            break;
        end
    end

    simData = beamSimData;
    simData.R  = REqu;
    simData.x  = xEqu;
    simData.xi = xi_k;


end

function equiEqu = implicitEquilibriumFunction(psi_k, simPars, discPars, Ba, Bc, useExp)
    % Function with the eom terms needed for the equilibrium configuration
    % (i.e., without the dynamic terms / only gravity and strain)
    arguments
        psi_k    (:,:) double

        % beamSimPars object containing the simulation parameters
        simPars (1,1) beamSimPars

        % Parameters of the discrete beam
        discPars (1,1) beamParamsDiscrete

        % Selection matrix for allowed discrete deformations (6, nAllwd)
        Ba   (6, :)     double

        % Selection matrix for constrained / constant discrete deformations
        % (6, 6-nAllwd)
        Bc   (6, :)     double

        % Used retraction map: 0 = Cayley, 1 = Exponential Map
        useExp (1,1) logical
    end

    % Nr. of nodes, segments and allowed deformation DoFs
    nSeg   = size(psi_k, 2);
    nNodes = nSeg + 1;


    %% Relative Kinematics
    % Compute full discrete deformations xi from allowed and constrained
    % deformations psiA and psiC
    xi_k = Ba * psi_k + Bc * Bc.' * simPars.xiRef;

    % Beam Jacobians and forward kinematics
    if useExp
        [~, ~,~, g_ij_k] = beamRelFwdKinExp(xi_k * discPars.l, simPars.g0(:,:,1));
        [R_k, ~, ~, JBeam] = beamRelKinReducedCompleteExp(xi_k, g_ij_k, simPars.g0(:,:,1), discPars.l, Ba);
    else
        [~, ~,~, g_ij_k] = beamRelFwdKin(xi_k * discPars.l, simPars.g0(:,:,1));
        [R_k, ~, ~, JBeam] = beamRelKinReducedComplete(xi_k, g_ij_k, simPars.g0(:,:,1), discPars.l, Ba);
    end


    %% Node Forces
    % i.e., EOM terms that are multiplied with the Jacobian

    % Prepare external forces
    f_node_k_b = zeros(6, nNodes);
    f_node_k_s = zeros(6, nNodes);
    if ~isempty(simPars.f_node_b) && all(size(simPars.f_node_b) == [6, nNodes])
        f_node_k_b = simPars.f_node_b;
    end
    if ~isempty(simPars.f_node_s) && all(size(simPars.f_node_s) == [6, nNodes])
        f_node_k_s = simPars.f_node_s;
    end

    for iN = 2:nNodes
        f_node_k_b(:, iN) = ...
            ...% External body-fixed forces
            - f_node_k_b(:, iN) ...
            ... % External spatial forces
            + [
            R_k(:,:,iN).' * -f_node_k_s(1:3,iN)
            R_k(:,:,iN).' * -f_node_k_s(4:6,iN)
            ] ...
            ...% Gravity
            + simPars.g * [
            discPars.m_a(iN) * cross( discPars.x_a(:,iN), R_k(:,:,iN).' * [0;0;1] )
            discPars.mNode(iN)*R_k(:,:,iN).' * [0;0;1] 
            ];
    end


    %% Segment Forces
    % I.e., terms that are not multiplied with the Jacobian

    % Ext. forces and stress
    f_seg_k = discPars.CgenRed * Ba.' * (xi_k - simPars.xiRef);


    %% Equilibrium equation (Force balance)

    % Residual
    equiEqu = f_seg_k(:);

    % Only evaluate for nodes > 1, since Jacobian of the first node is
    % always zero
    for iN = 1:nNodes
        equiEqu = equiEqu + JBeam(:, :, iN).' * f_node_k_b(:, iN);
    end
end

function [R, x, g, g_ij] = beamRelFwdKinExp(xi, g0)
    %% Compute Beam Forward Kinematics from Relative Deformations
    % i.e., compute absolute positions and rotations for each node
    %
    % Inputs:
    %  xi   Array of relative deformations with dimensions (6, nSeg)
    %
    %       Important: Must correspond to discrete *updates*, not discrete
    %       *gradients*! I.e., tau(xi) = g_a^{-1}*g_{a+1} must hold!
    %
    %  g0   SE3 configuration element of the first node (4x4 matrix)
    %
    % Outputs:
    %   R  Array of node rotation matrices with dimensions (3,3,nNodes)
    %   x  Array of node position vectors with dimensions (3,nNodes)
    %   g  Array of SE3 matrices with dimensions (4,4,nNodes) that describe
    %   the *absolute* configurations of the nodes w.r.t. g0
    %
    %   g_ij Array of SE3 matrices with dimensions (4,4,nSeg) that
    %   describe the *relative* configurations between the nodes, i.e.,
    %        g_j = g_i * g_ij
    %   holds; the definition is
    %        g_ij = tau( xi_ij)
    %   (where xi is the relative *update* as described above).
    %   Index i corresponds to the update described by xi with index i.


    % Nr. of nodes and segments
    nSeg   = size(xi, 2);
    nNodes = nSeg + 1;

    % Compute absolute configuration from relative deformations
    g    = zeros(4,4,nNodes);
    g_ij = zeros(4,4,nSeg);

    % Absolute configuration of the first node is given with g0
    g(:,:,1) = g0;

    % Compute relative update between first two nodes
    g_ij(:,:,1) = expSE3( xi(:, 1) );

    for iN = 2:nNodes
        if iN < nNodes
            g_ij(:,:,iN) = expSE3( xi(:, iN) );
        end
        g(:,:,iN) = g(:,:,iN-1) * g_ij(:,:,iN-1);
    end

    % Get arrays of rotation matrices and position vectors
    [R, x] = RxFromSE3Matrix(g);
end

function [R, x, g, JBeam] = beamRelKinReducedCompleteExp(xi, g_ij, g0, l, Ba)
    %% Compute Beam Forward Kinematics (with reduced deformation DoFs) from Relative Deformations
    % i.e., compute absolute positions and rotations for each node
    %
    % Inputs: See below.
    %
    % Outputs:
    %   R  Array of node rotation matrices with dimensions (3,3,nNodes)
    %   x  Array of node position vectors with dimensions (3,nNodes)
    %   g  Array of SE3 matrices with dimensions (4,4,nNodes) that describe
    %   the *absolute* configurations of the nodes w.r.t. g0
    %
    %   JBeam: Array of Beam Jacobians of all nodes with dimensions
    %          (6,6, nAllwd*nSeg, nNodes)
    %

    arguments
        % Array of discrete deformation gradients with dimensions (6, nSeg)
        %   Important: Must correspond to discrete *gradients*, not discrete
        %   *updates*! I.e., tau(xi*l) = g_a^{-1}*g_{a+1} must hold!
        xi (6,:)

        g_ij (4,4,:)

        %  SE3 configuration element of the first node (4x4 matrix)
        g0 (4,4)

        % Segment length
        l (1,1)

        % Selection matrix for allowed discrete deformations (6, nAllwd)
        Ba   (6, :)
    end

    % Nr. of nodes and segments
    nSeg   = size(xi, 2);
    nNodes = nSeg + 1;
    nAllwd = size(Ba, 2);

    g = zeros(4,4,nNodes);

    % Absolute configuration of the first node is given with g0
    g(:,:,1) = g0;

    % Array holding all beam Jacobians
    JBeam = zeros(6, nAllwd*nSeg, nNodes);%, class(xi));

    for iN = 2:nNodes

        % Compute the relative update of the segment *before* the current
        % node, which is needed for both the absolute configuration and the
        % beam Jacobian
        %g_ij = caySE3( xi(:, iN-1) * l );

        % Compute absolute configuration of the current node
        g(:,:,iN) = g(:,:,iN-1) * g_ij(:,:,iN-1);

        for ii = 1:iN
            if ii == (iN-1)
                JBeam(:, (ii*nAllwd)-(nAllwd-1) : ii*nAllwd, iN) = ...
                    + l ...
                    * expRTDSE3( -xi(:, iN-1) * l ) * Ba;
            elseif ii < (iN-1)
                JBeam(:, (ii*nAllwd)-(nAllwd-1) : ii*nAllwd, iN) = ...
                    + lAdSE3Inv( g_ij(:,:,iN-1) ) ...
                    * JBeam(:, (ii*nAllwd)-(nAllwd-1) : ii*nAllwd, iN-1);
            else
                % Case for iN > iCol: Column = zero
                % (Do nothing: has already been initialized to zero)
            end
        end

    end

    % Get arrays of rotation matrices and position vectors
    [R, x] = RxFromSE3Matrix(g);
end