%% Startup File: Add file paths and check important dependencies

% Absolute file path of this function
funPath = mfilename("fullpath");

% Get the directory of the function = repository root path
rootPath = fileparts(funPath);

% Add important paths
addpath(genpath(fullfile(rootPath, "build")));
addpath(genpath(fullfile(rootPath, "src")));
addpath(genpath(fullfile(rootPath, "scripts")));
addpath(genpath(fullfile(rootPath, "simulation-studies")));
addpath(genpath(fullfile(rootPath, "tests")));
