function [gRef, xiRef] = beamSimRefConf(nSeg, params)
    %% Compute Beam Reference Configuration
    % Configuration of the undeformed beam
    % Here: Straight beam
    % Computed from (constant) discrete deformations along the beam.

    disp('Computing Reference Configuration...')

    % Segment length
    l = params.L / nSeg;

    % Configuration of first node
    %RRef0 = eul2rotm([0, 0, -pi/2]);

    % Rotation matrix for 90° rotation about global y axis;
    % local z axis now aligned with global x axis
    RRef0 = [
        0  0 1
        0  1 0
        -1 0 0
        ];

    xRef0 = [0;0;0];
    g0 = SE3Matrix(RRef0, xRef0);

    % Vector of constant deformations used for all segments
    xiConst = [0; 0; 0; 0; 0; 1];

    % Array of segment deformations
    xiRef = repmat(xiConst, [1, nSeg]);

    % Forward kinematics
    [~, ~, gRef] = beamRelFwdKin(xiRef*l, g0);
end