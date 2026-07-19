function beamModels = defineSimStudyBeamModels(options)
    %% Define beam models used for simulation studies
    % To be used as a common base for other simulation scripts.

    arguments
        options.solverConfig (1,1) beamSolverConfig = beamSolverConfig()
    end

    % Global solver configuration:
    % Either use the one given as optional argument or use the default
    % values given in the class definition
    solverConfig = options.solverConfig;

    %% Define Models

    beamModels(5) = beamSimModel;

    % AbsKin LGVI (Simo-Reissner)
    beamModels(1).modelName        = "AbsKin SR General";
    beamModels(1).solverConfig     = solverConfig;
    beamModels(1).reducedParams.Ba; % Not needed
    beamModels(1).reducedParams.Bc;
    beamModels(1).funHandle        = @(beamPars, simPars, simModel) ...
        beamMdlAbsKinLGVI_general_mex(simPars, beamPars, simModel.solverConfig, zeros(6,0), zeros(6,0));

    % RelKin VarInt Broyden with all 6 DoF (Simo-Reissner)
    beamModels(2).modelName        = "RelKin SR Broyden";
    beamModels(2).solverConfig     = solverConfig;
    beamModels(2).reducedParams.Ba = eye(6);
    beamModels(2).reducedParams.Bc = zeros(6,0);
    beamModels(2).funHandle        = @(beamPars, simPars, simModel) ...
        beamMdlRelKinVarInt_Broyden_mex(simPars, beamPars, simModel.solverConfig, simModel.reducedParams.Ba, simModel.reducedParams.Bc);

    % RelKin VarInt Broyden with 3 DoF (Kirchhoff)
    beamModels(3).modelName        = "RelKin KH Broyden";
    beamModels(3).solverConfig     = solverConfig;
    beamModels(3).reducedParams.Ba = [ eye(3); zeros(3)];
    beamModels(3).reducedParams.Bc = [ zeros(3); eye(3)];
    beamModels(3).funHandle        = @(beamPars, simPars, simModel) ...
        beamMdlRelKinVarInt_Broyden_mex(simPars, beamPars, simModel.solverConfig, simModel.reducedParams.Ba, simModel.reducedParams.Bc);

    % RelKin VarInt Recursive with all 6 DoF (Simo-Reissner)
    beamModels(4).modelName        = "RelKin SR Recursive";
    beamModels(4).solverConfig     = solverConfig;
    beamModels(4).reducedParams.Ba = eye(6);
    beamModels(4).reducedParams.Bc = zeros(6,0);
    beamModels(4).funHandle        = @(beamPars, simPars, simModel) ...
        beamMdlRelKinVarInt_Recursive_mex(simPars, beamPars, simModel.solverConfig, simModel.reducedParams.Ba, simModel.reducedParams.Bc);

    % RelKin VarInt Recursive with 3 DoF (Kirchhoff)
    beamModels(5).modelName        = "RelKin KH Recursive";
    beamModels(5).solverConfig     = solverConfig;
    beamModels(5).reducedParams.Ba = [ eye(3); zeros(3)];
    beamModels(5).reducedParams.Bc = [ zeros(3); eye(3)];
    beamModels(5).funHandle        = @(beamPars, simPars, simModel) ...
        beamMdlRelKinVarInt_Recursive_mex(simPars, beamPars, simModel.solverConfig, simModel.reducedParams.Ba, simModel.reducedParams.Bc);

end