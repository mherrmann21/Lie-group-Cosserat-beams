function xDot = simResComputeVelocities(tout, R, x, eta) %#codegen
    %% Compute Velocities from simulation output data
    %i.e, translational velocities in the inertial frame

    arguments
        tout (:, 1)     double % Time vector
        R    (3,3,:,:)  double % Array of rotation matrices
        x    (3,:,:)    double % Array of position vectors
        eta  (6,:,:)    double % Array of R6/se3 velocity vectors
    end

    % Get nr. of nodes and segments
    nNodes = size(R, 3);

    %% Compute velocities in inertial frame
    xDot = zeros(3, nNodes, length(tout));
    for iStep = 1:length(tout)
        for iN = 1:nNodes
            xDot(:, iN, iStep) = R(:,:,iN,iStep) * eta(4:end,iN,iStep);
        end
    end
end
