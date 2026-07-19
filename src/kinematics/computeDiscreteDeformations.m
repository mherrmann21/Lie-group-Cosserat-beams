function xi = computeDiscreteDeformations(g)
    % Compute the discrete deformations xi that correspond to the given
    % beam configuration.
    %  Important: The computed xi values correspond to discrete *updates*, 
    %  not discrete *gradients*! I.e., tau(xi) = g_a^{-1}*g_{a+1} holds!
    %  If the discrete "gradients" are required, the result must be divided
    %  by l in the parent function.
    %
    % Inputs:
    %   g:  Array of SE3 node configuration matrices with dimension 
    %       (4, 4, nNodes)
    %   x:  2D array of node position vectors with dimension (3, nNodes)
    %
    % Outputs:
    %   xi: 2D array of discrete deformations with dimension (6, nSeg)

    % Get number of discrete segments
    nSeg = size(g, 3) - 1;

    xi = zeros(6, nSeg);
    for iN = 1:nSeg
        xi(:,iN) = cayInvSE3( g(:,:,iN) \ g(:,:,iN+1) );
    end
end