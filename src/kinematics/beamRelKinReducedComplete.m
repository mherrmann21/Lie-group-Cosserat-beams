function [R, x, g, JBeam] = beamRelKinReducedComplete(xi, g_ij, g0, l, Ba)
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
    %          (6, nAllwd*nSeg, nNodes)
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
                    * cayRTDSE3( -xi(:, iN-1) * l ) * Ba;
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
