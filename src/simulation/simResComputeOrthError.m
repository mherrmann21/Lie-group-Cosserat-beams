function orthErrorR = simResComputeOrthError(R)
    %% Compute Orthogonality Error for simulation output data

    arguments
        R (3,3,:,:) double % Array of rotation matrices
    end

    % Get nr. of nodes and segments
    nNodes = size(R, 3);

    %% Compute orthogonality error of rotation matrices
    orthErrorR = zeros( nNodes, size(R,4) );
    for iStep = 1:size(R,4)
        for iN = 1:nNodes
            orthErrorR(iN, iStep) = orthError( R(:,:,iN,iStep));
        end
    end
end
