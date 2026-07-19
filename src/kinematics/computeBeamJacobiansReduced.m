function [JBeam, g_xi] = computeBeamJacobiansReduced(xi, l, Ba)
    %% Compute Beam Forward Kinematics (with reduced deformation DoFs) from Relative Deformations
    % i.e., compute absolute positions and rotations for each node
    %
    % Inputs: See below.
    %
    % Outputs:
    %   JBeam: Array of Beam Jacobians of all nodes with dimensions 
    %          (6, nAllwd*nSeg, nNodes)
    %   g_xi:  SE3 update matrices corresponding to the discrete
    %   deformations xi

    arguments
        % Array of discrete deformation gradients with dimensions (6, nSeg)
        %   Important: Must correspond to discrete *gradients*, not discrete
        %   *updates*! I.e., tau(xi*l) = g_a^{-1}*g_{a+1} must hold!
        xi (6,:)

        % Segment length
        l (1,1)

        % Selection matrix for allowed discrete deformations (6, nAllwd)
        Ba   (6, :)
    end

    % Nr. of nodes and segments
    nSeg   = size(xi, 2);
    nNodes = nSeg + 1;
    nAllwd = size(Ba, 2);

    g_xi = zeros(4,4,nSeg);

    % Array holding all beam Jacobians
    JBeam = zeros(6, nAllwd*nSeg, nNodes);

    for iN = 2:nNodes
        % Compute the relative update of the segment *before* the current
        % node, which is needed for both the absolute configuration and the
        % beam Jacobian
        g_xi(:,:,iN-1) = caySE3( xi(:, iN-1) * l );

        for ii = 1:iN
            if ii == (iN-1)
                JBeam(:, (ii*nAllwd)-(nAllwd-1) : ii*nAllwd, iN) = ...
                    + l ...
                    * cayRTDSE3( -xi(:, iN-1) * l ) * Ba;
            elseif ii < (iN-1)
                JBeam(:, (ii*nAllwd)-(nAllwd-1) : ii*nAllwd, iN) = ...
                    + lAdSE3Inv( g_xi(:,:,iN-1) ) ...
                    * JBeam(:, (ii*nAllwd)-(nAllwd-1) : ii*nAllwd, iN-1);
            else
                % Case for iN > iCol: Column = zero
                % (Do nothing: has already been initialized to zero)
            end
        end

    end
end
