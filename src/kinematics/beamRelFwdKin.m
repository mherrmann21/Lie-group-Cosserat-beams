function [R, x, g, g_ij] = beamRelFwdKin(xi, g0)
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
    g_ij(:,:,1) = caySE3( xi(:, 1) );

    for iN = 2:nNodes
        if iN < nNodes
            g_ij(:,:,iN) = caySE3( xi(:, iN) );
        end
        g(:,:,iN) = g(:,:,iN-1) * g_ij(:,:,iN-1);
    end

    % Get arrays of rotation matrices and position vectors
    [R, x] = RxFromSE3Matrix(g);
end
