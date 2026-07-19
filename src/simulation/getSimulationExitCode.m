function exitCode = getSimulationExitCode(simData, metaDataSteps, simPars, solverConfig)
    %% Get Exit Code for Finished Beam Simulation
    %
    % This function generates an exit code that indicates whether the
    % simulation was successful:
    %
    % 0: Simulation failed (e.g., due to the solver not converging in some time
    %    steps)
    %
    % 1: Simulation successful, the specified solver tolerance was satisfied in
    %    all time steps / no tolerance violations
    %
    % 2: The max. nr. of implicit iterations was reached in at least one time
    %    step, but the residual is still within the error limit
    %    (beamSolverConfig.errorMarginLimit)
    %
    % Inputs: See arguments block below.
    %
    % Maximilian Herrmann
    % Chair of Automatic Control
    % TUM School of Engineering and Design
    % Technical University of Munich

    arguments
        simData         (1,1) beamSimData
        metaDataSteps   (1,1) beamSimMetaDataSteps
        simPars         (1,1) beamSimPars
        solverConfig    (1,1) beamSolverConfig
    end

    nSteps = round( simPars.tEnd / simPars.h );

    % Check if simulation was successful (=if any iterations failed)
    if length(simData.tout) < nSteps
        % Simulation failed
        exitCode = 0;
    else
        if ( ~any( metaDataSteps.ExitFlag(:) ) && ...
                max( metaDataSteps.ImplicitError(:)) < solverConfig.errorMarginLimit ...
                )
            % Everything fine, no tolerances violated
            exitCode = 1;
        elseif ( any( metaDataSteps.ExitFlag(:) ) && ...
                max( metaDataSteps.ImplicitError(:)) < solverConfig.errorMarginLimit ...
                )
            % Max. nr. of implicit iterations reached, but residual is
            % within error limit (beamSolverConfig.errorMarginLimit)
            exitCode = 2;
        else
            % Should not be reached
            error('Error computing simulation exit code.');
        end
    end
end
