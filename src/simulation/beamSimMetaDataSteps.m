classdef beamSimMetaDataSteps
    % beamSimMetaDataSteps class to store metadata for beam simulations
    % (for individual space/time steps)
    %
    % Maximilian Herrmann
    % Chair of Automatic Control
    % TUM School of Engineering and Design
    % Technical University of Munich
    %
    % Abbreviations used below:
    % nNodes:  Nr. of nodes (= nSeg + 1)
    % nSeg:    Nr of segments
    % nSteps:  Nr of time steps (of the simulation)

    properties
        % Orthogonality error of R, dimensions: (nNodes, nSteps)
        orthErrorR          (:,:) double

        % Nr of iterations of the implicit solver
        % Dimensions:
        % (nNodes, nSteps) for absKin Models
        % (1, nSteps) for RelKin Models
        ImplicitIterations  (:,:) double

        % Residual error of the DEL equations after convergence
        % Dimensions as ImplicitIterations
        ImplicitError       (:,:) double

        % Exit flag of the implicit solver (Note: Meaning can depend on the
        % used solver)
        % Dimensions as ImplicitIterations
        ExitFlag            (:,:) double
    end
end
