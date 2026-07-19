function E = computeBeamEnergyEvolution(x, eta, xi, xiRef, beamPars, simPars)
    %% Compute energy evolution of beam simulation results
    %
    % Calculates the evolution of the total energy H, which is composed of:
    % * Kinetic energy T (at each time step)
    % * Potential energy U (at each time step)
    % * Dissipated energy D (total energy up to the current time step)
    % ( * Input/external energy E (total energy up to the current time
    % step) )
    % from simulation results.
    %
    % The function uses absolute quantities (velocities and positions) to
    % evaluate the energy.
    %
    % Inputs:  See arguments block below.
    %
    % Outputs: Object with arrays of the computed energy evolutions as
    %          properties
    %
    % Maximilian Herrmann
    % Chair of Automatic Control
    % TUM School of Engineering and Design
    % Technical University of Munich

    arguments
        % Array of node positions with dimensions (3, nNodes, nSteps)
        x       (3,:,:) double

        % Array of node velocity vectors with dimensions (6, nNodes, nSteps)
        eta     (6,:,:) double

        % Array of segment deformations with dimensions (6, nSeg, nSteps)
        xi      (6,:,:) double

        % Array of reference deformations with dimensions (6, nSeg)
        xiRef   (6,:)   double

        % Standard beam parameter struct
        beamPars  (1,1)   beamParams

        simPars (1,1) beamSimPars
    end

    % Get number of discrete nodes and segments
    nNodes = size(eta, 2);
    nSeg = nNodes - 1;

    % Compute segment length
    l = beamPars.L / nSeg;


    %% Compute Energies
    % Compute all energies node/segment-wise for each time step;
    % use absolute / global kinematic quantities

    E = beamSimEnergies;
    E.T = zeros(size(x,3),1);
    E.U = zeros(size(x,3),1);
    E.V = zeros(size(x,3),1);
    E.H = zeros(size(x,3),1);

    fac = [0.5; ones(nSeg-1, 1); 0.5];

    for iStep = 1:size(x,3)
        for iN = 1:nNodes

            % Kinetic and Potential energy
            E.T(iStep) = E.T(iStep) ...
                + l / 2  ...
                * fac(iN) * eta(:,iN,iStep)' * beamPars.Mgen ...
                * eta(:,iN,iStep);

            E.U(iStep) = E.U(iStep) ...
                + l * fac(iN) ...
                * beamPars.m * simPars.g * [0;0;1]' * x(:,iN,iStep);

            % Strain energy
            if iN < nNodes
                E.V(iStep) = E.V(iStep) ...
                    + 1 / 2 * l ...
                    * ( xi(:,iN,iStep) - xiRef(:,iN) )' * beamPars.Cgen ...
                    * ( xi(:,iN,iStep) - xiRef(:,iN) );
            end
        end

    end

    % Normalize Potential energy: Set initial value to 0
    E.U = E.U - E.U(1);

    % Total energy
    E.H = E.T + E.U + E.V;
end
