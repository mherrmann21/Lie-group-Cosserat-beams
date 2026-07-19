%% Build all .mex functions for the beam simulation framework

clear


%% Code Generation Settings
% Store in configuration object of class 'coder.MexCodeConfig'.

cfg = coder.config('mex');
cfg.TargetLang = 'C++';
cfg.GenerateReport = false;

% Function inlining "always" seems improve runtime significantly
cfg.InlineBetweenMathWorksFunctions = 'Always';
cfg.InlineBetweenUserAndMathWorksFunctions = 'Always';
cfg.InlineBetweenUserFunctions = 'Always';

% Disable memory checks etc. for better runtime performance
cfg.IntegrityChecks = false;
cfg.ResponsivenessChecks = true;

% For profiling
%cfg.EnableMexProfiling = true;
%cfg.TargetLang = 'C';

% Build folder
targetDir = fullfile(getRootFolder, "build");


%% Build functions

fprintf("Compiling MEX functions...\n\n");

functionNames = [
    ... % Discrete Dynamic Models
    "beamMdlAbsKinLGVI_general"
    "beamMdlAbsKinLGVI_perNode"
    "beamMdlRelKinVarInt_Broyden"
    "beamMdlRelKinVarInt_Recursive"
    ... % Continuous-Time Model
    "beamMdlAbsKinCont_RHS"
    ... % Static Model
    "computeStaticEquilibriumRelKin"
    ... % Simulation / Evaluation Functions
    "computeBeamEnergyEvolution"
    "simResComputeDeformations"
    "simResComputeOrthError"
    "simResComputeVelocities"
    ];

for iFun = 1:numel(functionNames)
    fprintf("Compiling function ""%s""...\n", functionNames(iFun));
    codegen("-d", targetDir, "-o", ...
        fullfile(targetDir, functionNames(iFun) + "_mex"), ...
        "-config", cfg, functionNames(iFun));
end

%% End Script

disp('Code generation finished.')
