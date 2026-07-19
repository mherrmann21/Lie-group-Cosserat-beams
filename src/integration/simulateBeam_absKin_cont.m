function [simData, metaDataSteps] = simulateBeam_absKin_cont(beamPars, simPars, solverFunHandle, solverOptions)
    %% Simulate absKin continuous model
    % with given beam and simulation parameters

    arguments
        % Beam parameter struct
        beamPars         (1,1) beamParams

        % simPars object (struct) with simulation parameters.
        % See class definition for details.
        simPars         (1,1) beamSimPars

        % Function handle of the ODE solver to use (e.g., @ode15s)
        solverFunHandle (1,1) function_handle

        solverOptions   (1,1) struct = odeset()
    end


    %% Assign parameters to local variables

    nNodes = size(simPars.g0, 3);
    nSeg = nNodes - 1;

    tspan   = [0, simPars.tEnd];
    [R0,x0] = RxFromSE3Matrix(simPars.g0);

    % Initial condition state vector
    xState0 = [
        reshape(R0, [9*nNodes, 1]);
        reshape(x0, [3*nNodes, 1]);
        zeros(6*nNodes, 1)
        ];

    
    %% Integrate model

    % Show progress?
    % figure('Name', 'ODE integration progress'); % for odeplot
    % solverOptions = odeset( solverOptions, 'OutputFcn', @odeplot);

    % Function for the RHS of the function
    odeFun = @(t,y) beamMdlAbsKinCont_RHS_mex(y, simPars, beamPars);

    [tout, xout] = solverFunHandle( odeFun, tspan, xState0, solverOptions);


    %% Post-Process Results

    simData = beamSimData;
    simData.tout = tout;

    % Reassign output data to get proper dimensions

    % Rotation matrices: (3, 3, nNodes, nSteps)
    RVec = xout(:, 1:9*nNodes);
    simData.R = permute( ...
        reshape( RVec, [length(tout), 3, 3, nNodes]), ...
        [2,3,4,1] ...
        );

    % Position vectors: (3, nNodes, nSteps)
    xVec = xout(:, 9*nNodes+1: 12*nNodes);
    simData.x = permute( ...
        reshape( xVec, [length(tout), 3, nNodes]), ...
        [2,3,1] ...
        );

    % Velocity vectors: (6, nNodes, nSteps)
    etaVec = xout(:, end-6*nNodes+1:end);
    simData.eta = permute( ...
        reshape( etaVec, [length(tout), 6, nNodes]), ...
        [2,3,1] ...
        );

    % Compute deformations
    simData.xi = simResComputeDeformations_mex( ...
                simData.tout, simData.R, simData.x, simData.eta ...
                ) / (beamPars.L / nSeg);

    % Define step metadata object to be consistent with the general
    % function output structure
    metaDataSteps = beamSimMetaDataSteps;

    % Add implicit error value for compatibility with discrete integrator
    % functions
    metaDataSteps.ImplicitError = 0;

end

