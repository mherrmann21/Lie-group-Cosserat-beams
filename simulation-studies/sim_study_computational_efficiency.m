%% Simulation Time Study for Geometrically Exact Beam Models
% Execute beam simulations for various parameter settings and measure
% simulation times.
%
% * Find coarsest time step for convergence for all segment numbers
% * For one time step, compute simulation times for all segment numbers
%
% Maximilian Herrmann
% Chair of Automatic Control
% TUM School of Engineering and Design
% Technical University of Munich


%% Preparation

%#ok<*UNRCH>
clear
close all

% Disable warnings for nearly singular matrices; this can occur for badly
% conditioned cases (but is not important)
warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'Coder:MATLAB:illConditionedMatrix');
warning('off', 'MATLAB:illConditionedMatrix');


%% Script Settings

% Debug flag: Sets end time to small number and skips parameter loops
DEBUG_SCRIPT = 1;

% Plot and save the results of the individual simulations (for each model)?
PLOT_RESULTS_INDIV = 0;

% Save simulation data of the individual simulations (for each model)?
SAVE_DATA = 0;

% Generate and save animation of each material/geometry (from ref. model)
ANIMATE_RESULTS = 1;

% Directory where the results subfolder will be created
opts.resultsDir = fullfile(getRootFolder, "results", "runs");

% Suffix of the created subfolder (after date and time)
opts.folderSuffix = "_simResults_simTimes";


%%  Simulation Parameters

simPars = beamSimPars();

% End time / length
simPars.tEnd = 0.1;

% No gravity
simPars.g = 0;

% Global solver settings
solverConfig = beamSolverConfig;
solverConfig.errorMargin = 1e-11;
solverConfig.maxIterations = 10;
solverConfig.JacobianIterationThreshold = 5;
solverConfig.UseExactJacobian = true;

% External forces (at beam tip)
fMax = [0.5 0 0 0 2 2 ]' * 0.5;  % Max. force
fTEnd = 0.05;                    % Force impulse end time


%% Define models to simulate

beamModels = defineSimStudyBeamModels('solverConfig', solverConfig);


%% Simulation Case Parameters

% Step time
%h = 2.^ceil(linspace(-8, -19, 9));
%h = 2.^-18; % Test for time in space
%h = 2.^-[10:15, 15.5, 16, 16.5, 17, 18];
h = 2.^-(10:0.2:22);

% Nr. of segments n
nSeg = 2.^(2:6);

% Step time, at which all segment numbers are simulated
hAll = 2^-21;

% Beam materials and parameters
beamPars(1) = beamParams_mbsd_stiff_rod;
beamPars(1).d = zeros(6,1);

beamPars(2) = beamParams_mbsd_stiff_rod;
%beamPars(2).d = [ones(3,1)*1e-4; ones(3,1)*1e-2];
beamPars(2).d = ones(6,1)*1.2e-4;

% Reference simulation parameters
hRef    = 2^-22;
nSegRef = 128;

% Comparison Metric Parameters
nSegRefError = nSeg(1); %nan;
hRefError    = h(1); %nan;



%% Handle Debugging Flag

if DEBUG_SCRIPT

    % Set parameter vectors to scalar values
    beamPars = beamPars(1);
    h    = h(4);
    nSeg = nSeg(2);

    % Set short end time
    simPars.tEnd = 0.01;
end

%% Prepare output folder structure

[saveDirCase, saveDirAll] = prepareSimStudyFolderStructure(opts,beamPars,h,nSeg);


%% Turn on diary / write console output to file

diary(fullfile(saveDirAll, 'simulation.log'));

simStartTime = datetime;
fprintf('Starting Log. Time: %s\n', string(simStartTime, 'dd.MM.yy, HH:mm:ss'));

% Print host name to be able to identify PC afterwards
% Only works for windows!
% https://www.mathworks.com/matlabcentral/answers/398048-how-to-get-the-name-of-the-computer-under-matlab#answer_317798
fprintf('   Host Machine: %s\n\n', getenv('COMPUTERNAME'));

fprintf('Parameter Values:\n');
fprintf(  '   Mat:  %d values\n', numel(beamPars));
fprintf( ['   h:    ', repmat('%.3E ', 1, numel(h)), '\n'], h);
fprintf( ['   nSeg: ', repmat('%d ', 1, numel(nSeg)), '\n'], nSeg);
fprintf('\n')

fprintf('Solver Settings: errorMargin = %.1E, maxIterations = %d, JacobianIterationThreshold = %d\n\n', ...
    solverConfig.errorMargin, ...
    solverConfig.maxIterations, ...
    solverConfig.JacobianIterationThreshold);


%% Run All Simulation Cases

% Counter variable; just to display progress
caseCount = 0;

% Initialize results variables
% Note: Index order follows standard convention: Mat, h, nSeg
simIterationsAvg = nan(length(beamPars), length(h), length(nSeg), length(beamModels));
simIterationsMax = nan(length(beamPars), length(h), length(nSeg), length(beamModels));
simTimeTotal     = nan(length(beamPars), length(h), length(nSeg), length(beamModels));
simExitCode      = nan(length(beamPars), length(h), length(nSeg), length(beamModels));
simError         = nan(length(beamPars), length(h), length(nSeg), length(beamModels));


%% Simulate Case Models
for iMat = 1:length(beamPars)

    %% Run Reference Simulation

    beamSimRef = beamSimulation;
    beamSimRef.simPars = simPars;
    [beamSimRef.simPars.gRef, beamSimRef.simPars.xiRef] = beamSimRefConf(nSegRef, beamPars(iMat));
    beamSimRef.simPars.g0   = beamSimRef.simPars.gRef;
    beamSimRef.simPars.xi0  = beamSimRef.simPars.xiRef;
    beamSimRef.simPars.h    = hRef;

    % Add external forcing at the tip
    % (Smooth impulse (cosine) with specified max. force and end time)
    beamSimRef.simPars.force_scaling_mode = 1;
    beamSimRef.simPars.f_node_s           = [zeros(6,nSegRef),fMax];
    beamSimRef.simPars.force_tEnd         = fTEnd;

    beamSimRef.beamPars = beamPars(iMat);
    beamSimRef.simModel = beamModels(1); % Use AbsKin LGVI

    beamSimRef = beamSimRef.simulateModel("accurateTiming", false);

    % Plot results and save plots
    saveDirRefPlots = fullfile(saveDirAll, sprintf("Refsim_Mat_%d", iMat));
    if ~isfolder(saveDirRefPlots)
        mkdir(saveDirRefPlots)
    end
    disp('Generating Reference Simulation Plots...')
    beamSimRef.plotSimResults("hPlot", 2^-13, ...
        "savePlotsJPEG", ~DEBUG_SCRIPT, ...
        "savePlotsFig", ~DEBUG_SCRIPT, ...
        "saveDir", saveDirRefPlots ...
        );


    %% Run other simulations
    for iSeg = 1:length(nSeg)
        for ih = 1:length(h)
            %% Prepare Simulations

            fprintf('\n\nStarting Case %d/%d/%d: Mat=%d h=%.2E nSeg=%d \n', ...
                iMat, ih, iSeg, iMat, h(ih), nSeg(iSeg));

            % Assign Input Parameters for Current Simulation Case
            simPars.h = h(ih);

            caseCount = caseCount + 1;
            fprintf('     Overall: Case %d / %d\n\n', ...
                caseCount, length(beamPars)*length(h)*length(nSeg) )

            % Compute Reference and Initial Configuration
            [simPars.gRef, simPars.xiRef] = beamSimRefConf(nSeg(iSeg), beamPars(iMat));
            %[simPars.g0, simPars.xi0]  = beamSimInitialConf(nSeg(iSeg), beamPars(iMat));
            simPars.g0  = simPars.gRef;
            simPars.xi0 = simPars.xiRef;

            %% Simulate Models

            for iSim = 1:length(beamModels)
                close all
                % Check if we simulate the current time step:
                % * if it is the one where we do all simulations
                % * if the model has not converged yet

                if h(ih) == hAll || ~any(simExitCode(iMat, 1:ih, iSeg, iSim))

                    %% Simulate Model

                    fprintf('\nSimulating Model %d/%d...\n', iSim, length(beamModels))

                    % Create model object
                    beamSim = beamSimulation;
                    beamSim.simPars  = simPars;
                    beamSim.beamPars = beamPars(iMat);
                    beamSim.simModel = beamModels(iSim);

                    % Add external forcing at the tip
                    % (Smooth impulse (cosine) with specified max. force and end time)
                    beamSim.simPars.force_scaling_mode = 1;
                    beamSim.simPars.f_node_s           = [zeros(6,nSeg(iSeg)),fMax];
                    beamSim.simPars.force_tEnd         = fTEnd;

                    % Simulate model
                    beamSim = beamSim.simulateModel("accurateTiming", true);

                    if beamSim.simRes.metaDataSim.exitCode
                        % Compute simulation errors
                        [simError(iMat, ih, iSeg, iSim), ~] = simResComputeReferenceErrors( ...
                            hRefError, nSegRefError, beamPars(iMat).L, ...
                            beamSim.simRes.simData, beamSimRef.simRes.simData ...
                            );
                        if isnan(simError(iMat, ih, iSeg, iSim))
                            warning("Error while computing simulation error: NaN returned.")
                        end

                        % Store iteration count, simulation time etc.
                        simIterationsAvg(iMat, ih, iSeg, iSim) = beamSim.simRes.metaDataSim.ImplicitIterations.mean;
                        simIterationsMax(iMat, ih, iSeg, iSim) = beamSim.simRes.metaDataSim.ImplicitIterations.max;
                        simTimeTotal    (iMat, ih, iSeg, iSim) = beamSim.simRes.metaDataSim.TotalTime;
                        simExitCode     (iMat, ih, iSeg, iSim) = beamSim.simRes.metaDataSim.exitCode;
                    end


                    %% Plot results (if enabled)
                    % Only successful simulations to save some time
                    if PLOT_RESULTS_INDIV
                        if beamSim.simRes.metaDataSim.exitCode

                            disp('Generating Individual Plots...')

                            beamSim.plotSimResults( ...
                                "savePlotsJPEG", true, ...
                                "savePlotsFig", true, ...
                                "saveDir", saveDirCase(iMat, ih, iSeg), ...
                                "computeEnergy", true, ...
                                "hPlot", max([h(ih), 2^-13])...
                                );
                        end
                    end

                    %% Save data (with or without full simulation data)
                    disp('Saving Individual Data...');

                    % Remove the detailed simulation results data;
                    % only keeps lightweight aggregate simulation metadata
                    if ~SAVE_DATA
                        beamSim = beamSim.clearSimData();
                    end

                    % Save results
                    beamSim.saveObject( ...
                        fullfile( ...
                        saveDirCase(iMat, ih, iSeg), ...
                        strcat("beamMdl ", beamSim.simModel.modelName, ".mat") ) ...
                        );
                else
                    disp('    Skipping simulation because model has already converged.')
                end %
            end % for iSim
        end % for iSeg


        %% Animate Beam Results for Visual Check
        % Once for each Material/Geometry;
        % use the absKin reference simulation

        if ANIMATE_RESULTS
            disp('Animating Results...')

            movieSavePath = fullfile(saveDirAll, ...
                sprintf('animation_refMdl_case_%d', iMat)...
                );

            try
                animateBeamSimRes( ...
                    beamSim.simRes.simData, ...
                    beamPars(iMat), true, ...
                    movieSavePath ...
                    );
            catch ME
                warning(ME.identifier, ...
                    'Error during animation:\n %s', ME.message);
            end
        end
    end % for ih

    %% Save reference sim data (with or without full simulation data)
    disp('Saving Individual Data (reference simulation)...');

    % Remove the detailed simulation results data;
    % only keeps lightweight aggregate simulation metadata
    if ~SAVE_DATA
        beamSimRef = beamSimRef.clearSimData();
    end

    % Save results
    beamSimRef.saveObject( ...
        fullfile(saveDirAll, sprintf("beamMdlRef_Mat_%d.mat", iMat) ) ...
        );

end % for iMat

%% Save overall data

disp('Saving Overall Results Data...')

try
    save( ...
        fullfile(saveDirAll, "simTimeStudyResults.mat"), ...
        'simIterationsAvg', 'simIterationsMax', 'simTimeTotal', 'simExitCode', 'simError', ...
        'beamPars', 'nSeg', 'h', 'beamModels', 'saveDirCase',...
        '-mat', '-v7.3' ...
        );
catch ME
    warning(ME.identifier, 'Could not save data:\n %s', ME.message);
end


%% End script
simStopTime = datetime;

fprintf(...
    'Finished. Time: %s, Total duration: %s (hrs/min/s)\n', ...
    string(simStopTime, 'dd.MM.yy, HH:mm:ss'), ...
    string( duration(simStopTime-simStartTime, 'Format', 'hh:mm:ss') ) ...
    );

% Turn off diary
diary('off')
