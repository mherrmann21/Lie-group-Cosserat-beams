classdef beamSimMetaDataSim
    % beamSimMetaDataSim class to store metadata for beam simulations
    % (for the entire simulation)
    %
    % Contains purely scalar values.
    %
    % Maximilian Herrmann
    % Chair of Automatic Control
    % TUM School of Engineering and Design
    % Technical University of Munich

    properties
        ImplicitIterations = struct(...
            'min', [], ...
            'max', [], ...
            'mean', [] ...
            );

        ImplicitError = struct(...
            'min', [], ...
            'max', [], ...
            'mean', [] ...
            );

        StepTimeMean     (1,1) double
        TotalTime        (1,1) double

        TotalIterations  (1,1) double

        % For variational integrators: Exit code
        %  1: Everything fine, no tolerances violated
        %  2: Max. nr. of implicit iterations reached, but residual is
        %     within error limit (beamSolverConfig.errorMarginLimit)
        %     (simulation can still be considered successful)
        %  0: Max. nr. of implicit iterations reached, residual is above
        %     error limit. Simulation failed.
        exitCode        (1,1) double     
    end
end
