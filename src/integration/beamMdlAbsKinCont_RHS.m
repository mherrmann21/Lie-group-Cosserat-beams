function F = beamMdlAbsKinCont_RHS(xState, simPars, params)
    %% Right-Hand Side of Cont.-Time Beam Model in Absolute Description
    %
    % Beam Model:
    % * Continuous-time
    % * Absolute kinematics formulation using rotation matrices to describe
    %   the configuration
    % * "Node-level" model with individual EOMs for each Node
    %
    % This function computes the right-hand side of the beam EOMs.
    %
    % Inputs: See arguments block below.
    %
    % Maximilian Herrmann
    % Chair of Automatic Control
    % TUM School of Engineering and Design
    % Technical University of Munich

    arguments
        % State vector / array of model state variables with length
        % 18 * (n+1) in the form [ vec(R), x, eta ]
        % (with se3 velocities eta as R6 vectors; n: nr. of segments)
        xState (:,1) double

        % Simulation parameters
        simPars  (1,1) beamSimPars

        % Beam parameter object
        params   (1,1) beamParams
    end


    % Get number of discrete nodes and segments
    nNodes = length(xState) / 18;
    nSeg = nNodes - 1;

    % Segment length
    l = params.L / nSeg;

    % Get configuration arrays from state vector
    [RVec, xVec, etaVec] = getStateVectorComponents(xState);
    R   = reshape(RVec, [3,3,nNodes]);
    x   = reshape(xVec, [3, nNodes]);
    eta = reshape(etaVec, [6, nNodes]);
    g   = SE3Matrix(R,x);

    % Compute discrete deformations
    % Note: The computed values for xi correspond to discrete updates, not
    % gradients; thus, they are divided by segment length l afterwards
    xi = computeDiscreteDeformations(g) / (params.L / nSeg);

    % Compute strain rates (for dissipation)
    xi_dot = zeros(6, nSeg);
    for iN = 1:nSeg
        xi_dot(:,iN) = cayRTDInvSE3(-l*xi(:,iN)) * (...
            eta(:,iN+1) - lAdSE3Inv(caySE3(l*xi(:,iN)))*eta(:,iN)...
            ) / l;
    end


    %% Boundary conditions: Set pose of first node to zero (cantilever beam)
    isCantilever = 1;%~all(isnan(simPars.g0));
    if isCantilever
        % Fix configuration of first node to given values
        R(:,:,1) = simPars.g0(1:3,1:3,1);
        x(:,1)   = simPars.g0(1:3,4,1);
    end


    %% Compute RHS for all nodes

    Fvel = zeros(6, nNodes);

    % If the beam is a cantilever: Skip first node and set its velocities
    % to zero
    if isCantilever
        nStart = 2;
        Fvel(:, 1) = zeros(6,1);
    else
        nStart = 1;
    end

    % Pre-compute stress and damping force for all nodes (in the "relative"
    % space)
    tauStrDmp = ...
        + params.Cgen * (xi - simPars.xiRef) ...
        + diag(params.d)*xi_dot;

    for iN = nStart:nNodes

        %%% Strain term: Stiffness and Damping
        if iN == 1
            % First Node
            fStrDmp_i = -2 / l * cayRTDInvSE3(xi(:, iN) * l).' * tauStrDmp(:, iN);
        elseif iN == nNodes
            % Last Node
            fStrDmp_i = + 2 / l * lAdSE3( caySE3(xi(:, iN-1) * l) ).' * ...
                cayRTDInvSE3(xi(:, iN-1) * l).' * tauStrDmp(:, iN-1);
        else
            % Interior nodes
            fStrDmp_i = ...
                - 1/l * cayRTDInvSE3(xi(:, iN) * l).' * tauStrDmp(:, iN) ...
                + 1/l * lAdSE3( caySE3( xi(:, iN-1) * l ) ).' ...
                * cayRTDInvSE3(xi(:, iN-1) * l).' * tauStrDmp(:, iN-1);
        end

        %%% Compute complete right-hand side
        Fvel(:, iN) = params.MgenInv * ( ...
            + CoSadSE3( eta(:,iN) ) * params.Mgen * eta(:,iN) ...
            - [ zeros(3,1); R(:,:,iN).' * params.m * simPars.g * [0;0;1]] ...
            - fStrDmp_i ...
            );

        % For end node: Add constant external forces
        % (Both body-fixed and spatial)
        % if iN == nNodes
        %     Fvel(:, iN) =  Fvel(:, iN) + params.MgenInv * (...
        %         + 2/l* simPars.f_tip_b ...
        %         + 2/l* [ ...
        %         R(:,:,end).' * simPars.f_tip_s(1:3);
        %         R(:,:,end).' * simPars.f_tip_s(4:6)
        %         ] ...
        %         );
        % end
    end


    %%% Compute RHS of configurations / kinematic equations
    [RDot, xDot] = beamAbsKinematicsRHS(R, eta);

    % Boundary conditions: Set velocities of first node to zero (cantilever beam)
    if isCantilever
        RDot(:,:,1) = zeros(3,3);
        xDot(:,1)   = zeros(3,1);
    end

    %%% Assemble full RHS equation
    F = [
        reshape(RDot, [9*nNodes, 1]);
        reshape(xDot, [3*nNodes, 1]);
        reshape(Fvel, [6*nNodes, 1])
        ];

end


%% Local functions

function [RDot, xDot] = beamAbsKinematicsRHS(R,eta)
    % Compute the right-hand side of the beam kinematics equation; i.e.,
    % the kinematics equations for SO3 and R3
    %
    % Inputs:
    %   R:  3D array of node rotation matrices with dimension (3, 3, nNodes)
    %  eta: 2D array of node velocities with dimension (6, nNodes)
    %
    % Output:
    %  RDot: 3D array of node rotation matrix derivatives with dimension
    %        (3, 3, nNodes)
    %  xDot: 2D array of node position derivative vectors with dimension
    %        (3, nNodes)

    % Get number of discrete nodes ( = Nr. of segments +1)
    nNodes = size(R, 3);

    RDot = zeros(3,3, nNodes);
    xDot = zeros(3, nNodes);
    for iN = 1:nNodes
        % Rotational and translational kinematic equations
        RDot(:,:,iN) = R(:,:,iN) * skew( eta(1:3, iN) );
        xDot(:,iN)   = R(:,:,iN) * eta(4:end, iN);
    end
end
