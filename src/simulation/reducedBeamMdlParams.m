classdef reducedBeamMdlParams
    % Class to store parameters (selection matrices) for reduced relative
    % beam models
    % Note: default values define a Kirchhoff beam

    properties
        Ba (:,:) double = [eye(3); zeros(3)];
        Bc (:,:) double = [zeros(3); eye(3)];
    end
end
