function [g0, xi0] = beamSimInitialConf(nSeg, params)
    %% Compute Initial Beam Configuration
    
    disp('Computing Initial Configuration...')

    % Segment length
    l = params.L / nSeg;

    % Configuration of first node
    R00 = eul2rotm([0, 0, -pi/2]);
    x00 = [0;0;0];
    g00 = SE3Matrix(R00, x00);

    % Vector of constant deformations used for all segments
    xiRot = [pi/2; pi/2; 0.2];
    xiTrl = [0; 0; 1];
    xiConst = [ xiRot; xiTrl ];

    % Array of segment deformations
    xi0 = repmat(xiConst, [1, nSeg]);%.* repmat(linspace(1,0,nSeg), [6,1]).^2;
    xi0(end,:) = 1;


    %xi0(:,end) = [zeros(5,1);1];

    % Forward kinematics
    [~, ~, g0] = beamRelFwdKin(xi0*l, g00);
end
