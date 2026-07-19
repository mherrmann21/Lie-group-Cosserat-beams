classdef beamSolverConfig
    % Class to hold the (implicit) solver configuration for a beam model

    properties
        % Target value for the solver error margin
        errorMargin                 (1,1) double = 1e-8;

        % Solver error margin at which the simulation is cancelled
        errorMarginLimit            (1,1) double = 1e-8;

        maxIterations               (1,1) double = 100;

        % For relKin Broyden and absKin general integrator:
        % Nr. of iterations that are allowed in one time step before the
        % Jacobian matrix is recomputed
        JacobianIterationThreshold  (1,1) double = 4;

        % For absKin general integrator:
        % In the dissipative case, use the exact Jacobian (containing the
        % additional dissipation term; slow to compute) or only the
        % Jacobian from the non-dissipative case (faster to compute)
        % (Exact Jacobian necessary if strong dissipation is present)
        UseExactJabocobian          (1,1) logical = false;
    end
end
