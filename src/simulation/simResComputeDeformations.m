function xi = simResComputeDeformations(tout, R, x, eta)
    %% Compute discrete deformations for simulation output data
    %
    %  Important: The computed xi values correspond to discrete *updates*, 
    %  not discrete *gradients*! I.e., tau(xi) = g_a^{-1}*g_{a+1} holds!
    %  If the discrete "gradients" are required, the result must be divided
    %  by l in the parent function.

    arguments
        tout (:, 1)     double % Time vector
        R    (3,3,:,:)  double % Array of rotation matrices
        x    (3,:,:)    double % Array of position vectors
        eta  (6,:,:)    double % Array of R6/se3 velocity vectors
    end

    % Get nr. of nodes and segments
    nNodes = size(R, 3);
    nSeg   = nNodes - 1;


    %% Compute discrete deformations
    xi = zeros(6, nSeg, length(tout));
    for iStep = 1:length(tout)
        xi(:,:,iStep) = computeDiscreteDeformations(SE3Matrix(R(:,:,:,iStep), x(:,:,iStep) ));
    end

end
