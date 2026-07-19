%% Add repository folders to the MATLAB path

% Absolute path of this startup script
scriptPath = mfilename("fullpath");

% The startup script is located in the repository root
rootPath = fileparts(scriptPath);

% Add important paths
addpath(genpath(fullfile(rootPath, "build")));
addpath(genpath(fullfile(rootPath, "src")));
addpath(genpath(fullfile(rootPath, "scripts")));
addpath(genpath(fullfile(rootPath, "simulation-studies")));
addpath(genpath(fullfile(rootPath, "tests")));
