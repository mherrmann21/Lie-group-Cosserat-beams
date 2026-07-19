classdef beamSimModel
    %beamSimModel
    % Class to store all relevant properties of a given beam model,
    % including solver settings etc.

    properties
        % Function handle to the simulateBeam_ function for the given model
        funHandle       (1,1) function_handle = @beamMdlRelKinVarInt_Broyden_mex;

        % Name for the model
        modelName       (1,1) string = "RelKin Broyden";

        solverConfig    (1,1) beamSolverConfig

        reducedParams   (1,1) reducedBeamMdlParams
    end
end
