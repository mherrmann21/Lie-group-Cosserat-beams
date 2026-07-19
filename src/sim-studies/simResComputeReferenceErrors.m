function [errorRel, errorMat] = simResComputeReferenceErrors(hEval, nSegEval, L, simDataComp, simDataRef)
    %% Compute Error Metric of a Beam Simulation Compared to a Reference Simulation
    %
    % Method after [Dem+15, pp.113]; difference is that here we used the
    % configuration difference between to space-time nodes expressed with
    % the relative se3 element between them as the error measure.
    %
    %  * The comparison is done at the space-time-grid defined by
    %    hEval, nSegEval
    %  * The data of reference and comparison simulation is resampled to
    %    that grid by simple LINEAR interpolation
    %    -> Not good for rotation matrices
    %    -> ToDo: Implement proper Lie group interpolation on SE(3)
    %
    %
    % Inputs: See below.
    %
    % Outputs:
    %   * errorRel: Scalar relative error of the overall simulation
    %   * errorMat: Error matrix of all space-time nodes
    %
    % Maximilian Herrmann
    % Chair of Automatic Control
    % TUM School of Engineering and Design
    % Technical University of Munich

    arguments
        % Time step and nr. of segments that define the space-time grid,
        % at which the reference error is evaluated
        hEval       (1,1) double
        nSegEval    (1,1) double

        % Beam length
        L           (1,1) double

        % simData object (with simulation results data) of the comparison
        % and the reference simulation
        % (needed fields: .R, .x, .tout)
        simDataComp  (1,1) beamSimData
        simDataRef   (1,1) beamSimData
    end


    %% Interpolate data to required grid
    fprintf('Computing Simulation Error Metrics...\n');
    fprintf('   Used space-time grid: h=%.2E, nSeg=%d\n', hEval, int32(nSegEval));

    tEnd = min([simDataComp.tout(end), simDataRef.tout(end)]);

    lVecComp = linspace(0, 1, size(simDataComp.R, 3));
    lVecRef  = linspace(0, 1, size(simDataRef.R, 3));
    lVecEval = linspace(0, 1,nSegEval + 1);

    % Time vector of the evaluation (query) grid
    tVecEval = 0:hEval:tEnd;


    % Comparison simulation
    RComp = interpn( ...
        (1:3).', (1:3).', lVecComp, simDataComp.tout, ...
        simDataComp.R, ...
        (1:3).', (1:3).', lVecEval, tVecEval);

    xComp = interpn( ...
        (1:3)', lVecComp, simDataComp.tout, ...
        simDataComp.x, ...
        (1:3)', lVecEval, tVecEval);

    % Reference simulation
    RRef = interpn( ...
        (1:3)', (1:3)', lVecRef, simDataRef.tout, ...
        simDataRef.R, ...
        (1:3)', (1:3)', lVecEval, tVecEval');

    xRef = interpn( ...
        (1:3)', lVecRef, simDataRef.tout, ...
        simDataRef.x, ...
        (1:3)', lVecEval, tVecEval);


    %% Compute relative error
    % Note: Exclude first and last time step

    % 0 = based on retraction map, 1 = as in [Dem+15]
    ERROR_TYPE = 0;

    % Initialize error matrix of all space-time nodes
    errorMat = nan(length(tVecEval)-2, nSegEval+1);

    for iStep = 2:(length(tVecEval)-1)
        for iN = 1:(nSegEval+1)

            switch ERROR_TYPE
                case 0
                    % Get SE3 matrices for the current node for both simulations
                    gRefInv = SE3MatrixInv( ...
                        RRef(:,:, iN, iStep), ...
                        xRef(:,   iN, iStep) ...
                        );

                    gSim = SE3Matrix( ...
                        RComp(:,:, iN, iStep), ...
                        xComp(:,   iN, iStep) ...
                        );

                    % Compute se3 element corresponding to the relative configuration
                    % between them
                    xi_err = cayInvSE3( gRefInv*gSim );

                    % Compute norm as error metric
                    errorMat(iStep-1, iN) = norm(xi_err);
                case 1
                    % Get SE3 matrices for the current node for both simulations
                    gRef = SE3Matrix( ...
                        RRef(:,:, iN, iStep), ...
                        xRef(:,   iN, iStep) ...
                        );

                    gSim = SE3Matrix( ...
                        RComp(:,:, iN, iStep), ...
                        xComp(:,   iN, iStep) ...
                        );

                    % Compute norm as error metric
                    errorMat(iStep-1, iN) = ...
                        norm( gSim-gRef) / norm(gRef);
                otherwise
            end

        end
    end

    % Compute relative error of the entire Simulation via Frobenius norm of the
    % error matrix (= MATLAB 'norm' for matrices)
    % Formula from [Dem+15, p.115]
    errorRel = norm(errorMat) * sqrt(hEval) * sqrt(L/nSegEval);

    fprintf('   Relative Error: %e\n', errorRel);
end

