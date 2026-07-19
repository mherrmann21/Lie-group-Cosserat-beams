classdef beamSimPars
    % beamSimPars class containing all data for a specific simulation
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
        g0      (4,4,:) double   % SE3 node configuration matrices for the initial beam conf. (4,4,nNodes)
        xi0     (6,:)   double   % Discrete deformations for the initial beam conf.           (6, nSeg)

        gRef    (4,4,:) double   % Same as g0, but for the beam reference conf.
        xiRef   (6,:)   double   % ...

        tEnd    (1,1)   double   % simulation end time (= length of the simulation)
        h       (1,1)   double   % Sampling time (for discrete models) (s)

        g       (1,1)   double  = 9.81;  % Gravity constant


        %% Time-dependent External forces

        % External node forces (wrenches) in body-fixed (b) and spatial (s) frames
        f_node_b      (6,:) double
        f_node_s      (6,:) double

        % Time scaling for external forces
        % 0 = Constant
        % 1 = Smooth impulse that begins at t=0 and ends at t = tEnd
        % 2 = Smooth increase from 0 at t=0 to max. force at t = tEnd
        % 3 = Smooth decrease from max force at t=0 to 0 at t = tEnd
        % (4 = arbitrary linear interpolation with given time/scaling
        % vectors)
        force_scaling_mode    (1,1) int8

        % End time of transient part
        force_tEnd      (1,1) double

        % Vector of time sample points for the time-dependent scaling
        f_node_tVec     (1,:) double
        
        % Vector of scale values (at the time sample points)
        f_node_sVec     (1,:) double
        

        %% Additional rigid bodies attached to the beam nodes

        % Array of generalized inertia tensors (w.r.t. the beam center 
        % line/body-fixed cross-section frame), dimensions (6,6,nNodes)
        M_a         (6,6,:) double

        % Vector of body masses; dimensions (nNodes,1)
        m_a         (:,1)   double

        % Array of position vectors from the cross-section frame to the COM
        % frame of attached rigid bodies
        x_a         (3,:) double
    end
end
