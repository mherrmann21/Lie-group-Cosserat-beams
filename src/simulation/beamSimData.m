classdef beamSimData
    % simData class to store simulation data for beam simulations
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
        eta  (6,:,:)   double  % se(3)/ R6 node velocity vectors         (6, nNodes, nSteps)
        mu   (6,:,:)   double  % discrete momentum                       (6, nNodes, nSteps)
        R    (3,3,:,:) double  % node rotation matrices                  (3, 3, nNodes, nSteps)
        x    (3,:,:)   double  % node position vectors                   (3, nNodes, nSteps)
        xi   (6,:,:)   double  % segment deformations                    (6, nSeg, nSteps)
        xDot (3,:,:)   double  % node velocities in the inertial frame   (3, nNodes, nSteps))
        tout (:,1)     double  % time values                             (nSteps, 1)
    end
end
