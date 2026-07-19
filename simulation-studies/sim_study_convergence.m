%% Convergence Study for Geometrically Exact Beam Models
% Compare convergence in :
%   * space discretization (for fixed time step)
%  or
%   * time discretization (for fixed space step)
%
% Methodology follows [Dem+15].
% We compute the simulation error metric for two reference simulations; one
% for the absKin LGVI and one for the relKinRed varInt model. This allows
% to both prove "pure" convergence of all models (using a reference
% simulation of the same model type) as well as a comparison, of how close
% the reduced model is to the full model.
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

% Run convergence study in space or time? ("space" or "time"?)
%CONV_STUDY_TYPE = "space";
CONV_STUDY_TYPE = "time";

% Plot the results of the individual simulations (for each model)?
PLOT_RESULTS_INDIV = 1;

% Save simulation data of the individual simulations (for each model)?
SAVE_DATA = 0;

% Generate and save animation of each material/geometry (from ref. model)
ANIMATE_RESULTS = 1;

% Directory where the results subfolder will be created
opts.resultsDir = fullfile(getRootFolder, "results", "runs");

% Suffixes of the created subfolder (after date and time)
RESULTS_FOLDER_SUFFIX_SPACE = '_simResults_convSpace';
RESULTS_FOLDER_SUFFIX_TIME  = '_simResults_convTime';

%% Constant Simulation Parameters

simPars = beamSimPars();

% End time / length
simPars.tEnd = 0.1;

% No gravity
simPars.g = 0;

% External forces (at beam tip)
fMax = [0.5 0 0 0 2 2 ]'*10;    % Max. force
fTEnd = 0.05;                   % Force impulse end time

% Global solver settings (for relative and absolute models)
solverConfig = beamSolverConfig;
solverConfig.errorMargin = 1e-13;
solverConfig.maxIterations = 10;
solverConfig.JacobianIterationThreshold = 5;


%% Define models to simulate

beamModels = defineSimStudyBeamModels('solverConfig', solverConfig);

if DEBUG_SCRIPT
    beamModels = beamModels(5);
end

%% Define reference simulation models

refBeamModels(2) = beamSimModel;

% AbsKin LGVI (Simo-Reissner)
refBeamModels(1).modelName        = "Ref AbsKin SR";
refBeamModels(1).solverConfig     = solverConfig;
refBeamModels(1).reducedParams;   % Not needed
refBeamModels(1).funHandle        = @(beamPars, simPars, simModel) ...
    beamMdlAbsKinLGVI_general_mex(simPars, beamPars, simModel.solverConfig, zeros(6,0), zeros(6,0));

% RelKin VarInt with 3 DoF (Kirchhoff)
refBeamModels(2).modelName        = "Ref RelKin KH";
refBeamModels(2).solverConfig     = solverConfig;
refBeamModels(2).reducedParams.Ba = [ eye(3); zeros(3)];
refBeamModels(2).reducedParams.Bc = [ zeros(3); eye(3)];
refBeamModels(2).funHandle        = @(beamPars, simPars, simModel) ...
    ...%beamMdlRelKinReducedVarIntRecursive_mex(simPars, beamPars, simModel.solverConfig, simModel.reducedParams.Ba, simModel.reducedParams.Bc);
    beamMdlRelKinVarInt_Broyden_mex(simPars, beamPars, simModel.solverConfig, simModel.reducedParams.Ba, simModel.reducedParams.Bc);

if DEBUG_SCRIPT
    refBeamModels = refBeamModels(2);
end

%% Simulation Case Parameters
switch CONV_STUDY_TYPE
    case "space"
        %h= 2^-21; % Above CFL limit for 2^7 segments
        h= 2^-19; % Good enough
        nSeg = 2.^(2:6);

    case "time"
        h = 2.^-(16:20);
        nSeg = 2^5;
    otherwise
end

% Beam materials and parameters
beamPars(1) = beamParams_mbsd_soft_rod;
beamPars(1).d = zeros(6,1);

%beamPars(2) = beamParams_LLA11_steelString;
%beamPars(2).d = [ones(3,1)*1e-3; ones(3,1)*1];


%%% Parameters for the reference simulation ("ground truth")

% For space convergence
% (step time h of the reference simulation is always h of the current case)
nSegRef = 2^7;

% For time convergence
% (nSeg of the reference simulation is always nSeg of the current case)
hRef = 2^-21;


%% Comparison Metric Parameters
% Space-Time values at which to evalue the comparison metric.
% If set to NaN: The space-time grid of the current case is used (i.e.,
% different grid for each case)

nSegRefError = nSeg(1); %nan;
hRefError    = h(1); %nan;


%% Handle Debugging Flag

if DEBUG_SCRIPT
    % Short parameter vectors
    beamPars = beamPars(1);
    %h    = h(round(length(h)/2));
    %nSeg = nSeg(round(length(nSeg)/2));

    h = h(round(linspace(1,length(h),3)));
    % h = h(1);

    % Set short end time
    %simPars.tEnd = 0.05;
end

%% Prepare output folder structure

disp('Preparing Output Folder Structure...')

% Get subfolder for all simulation results
switch CONV_STUDY_TYPE
    case "space"
        opts.folderSuffix = RESULTS_FOLDER_SUFFIX_SPACE;
    case "time"
        opts.folderSuffix = RESULTS_FOLDER_SUFFIX_TIME;
    otherwise
        error("Error: Wrong convergence study type (Allowed: ""space"" or ""time"").");
end

[saveDirCase,saveDirAll] = prepareSimStudyFolderStructure(opts,beamPars,h,nSeg);


%% Turn on diary / write console output to file

diary(fullfile(saveDirAll, 'simulation.log'));

simStartTime = datetime;
fprintf('Starting Log. Time: %s\n', string(simStartTime, 'dd.MM.yy, HH:mm:ss'));
fprintf('Convergence study type: %s\n\n', CONV_STUDY_TYPE);

% Print host name to be able to identify PC afterwards
% Only works for windows, see
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

% Initialize simulation objects for reference simulations
beamSimRef(length(refBeamModels),1) = beamSimulation;

% Initialize results variables
% Note: Index order follows standard convention: Mat, h, nSeg
simErrors          = nan(length(beamPars), length(h), length(nSeg), length(beamModels), length(refBeamModels));
simTimeTotal       = nan(length(beamPars), length(h), length(nSeg), length(beamModels));


% Initialize variables for the outer and loop inner loop, depending on the
% convergence study type
% (The outer loop has the fixed parameter, i.e., the one that is kept
% constant between reference and model simulation)
switch CONV_STUDY_TYPE
    case "space"
        parOut = h;
        parIn  = nSeg;
    case "time"
        parOut = nSeg;
        parIn  = h;
    otherwise
        error("Error: Wrong convergence study type (Allowed: ""space"" or ""time"").");
end

%% Simulate Case Models
for iMat = 1:length(beamPars)
    for iOut = 1:length(parOut)
        %% Run Reference Simulations
        % where the analyzed parameter is kept constant; the other one is
        % varied

        switch CONV_STUDY_TYPE
            case "space"
                hRefCur    = h(iOut);
                nSegRefCur = nSegRef;
            case "time"
                hRefCur    = hRef;
                nSegRefCur = nSeg(iOut);
        end

        fprintf('\n\nStarting Reference Simulations for Mat=%d, h=%.2E=2^%.2f, nSeg=%d...\n\n', ...
            iMat, hRefCur, log2(hRefCur), nSegRefCur ...
            );

        checkCFLLimit(beamPars(iMat), nSegRefCur, 'h', hRefCur);

        simParsRef   = beamSimPars();
        simParsRef.h = hRefCur;
        simParsRef.tEnd = simPars.tEnd;

        simParsRef.g = 0;

        % Add external forcing at the tip
        % (Smooth impulse (cosine) with specified max. force and end time)
        simParsRef.force_scaling_mode = 1;     
        simParsRef.f_node_s           = [zeros(6,nSegRefCur),fMax];  
        simParsRef.force_tEnd         = fTEnd;

        % Compute Reference and Initial Configuration
        [simParsRef.gRef, simParsRef.xiRef] = beamSimRefConf(nSegRefCur, beamPars(iMat));
        %[simParsRef.g0, simParsRef.xi0] = beamSimInitialConf(nSegRefCur, beamPars(iMat));
        simParsRef.g0  = simParsRef.gRef;
        simParsRef.xi0 = simParsRef.xiRef;

        % Visualize initial and reference configuration
        visualizeBeamConfig( simParsRef.g0, beamPars(iMat), ...
            'Initial Configuration (Ref. Simulation)');

        visualizeBeamConfig( simParsRef.gRef, beamPars(iMat), ...
            'Reference Configuration (Ref. Simulation)');

        %%% Run reference simulations
        % Variable to indiciate whether any of the reference simulations
        % were successful
        anyRefSimSuccessful = false;

        for iSim = 1:length(refBeamModels)
            beamSimRef(iSim).simPars  = simParsRef;
            beamSimRef(iSim).beamPars = beamPars(iMat);
            beamSimRef(iSim).simModel = refBeamModels(iSim);

            % Simulate model
            beamSimRef = beamSimRef(iSim).simulateModel;

            % Update success variable
            anyRefSimSuccessful = ...
                anyRefSimSuccessful || beamSimRef(iSim).simRes.metaDataSim.exitCode;


            %% Debug Plot Solver Stats
            if DEBUG_SCRIPT
                nameStr = 'Ref Sim ';
                if ~isempty(beamSimRef(iSim).simRes.metaDataSteps.ImplicitIterations)
                    tout = beamSimRef(iSim).simRes.simData.tout;

                    figHandles.SolverStats = figure( ...
                        'Name', [nameStr, 'Solver Stats'], 'NumberTitle','off');

                    t = tiledlayout(figHandles.SolverStats, 2,2);

                    % Residual Error
                    ax = nexttile(t);
                    plot(ax, tout, beamSimRef(iSim).simRes.metaDataSteps.ImplicitError');
                    title(ax, 'Residual Error', 'interpreter', 'latex')
                    xlim(ax, [0, tout(end)])
                    grid on

                    % Implicit solver error flag
                    ax = nexttile(t);
                    plot(ax, tout, beamSimRef(iSim).simRes.metaDataSteps.ExitFlag', '-o', 'MarkerSize', 3)
                    title(ax, 'Implicit Solver Error Flag', 'interpreter', 'latex')
                    xlim(ax, [0, tout(end)])
                    grid on

                    % Nr. of iterations (of the implicit solver)
                    ax = nexttile(t, [1,2]);
                    plot(ax, tout, beamSimRef(iSim).simRes.metaDataSteps.ImplicitIterations', '-o');
                    title(ax, 'Nr. of Iterations', 'interpreter', 'latex')
                    xlim(ax, [0, tout(end)])
                    grid on

                    xlabel(ax, 'time $t$ / s', 'interpreter', 'latex')
                end
            end

        end

        %%% Check if any of the reference simulations were successful;
        % If not, we don't even need to run the individual models
        if ~anyRefSimSuccessful
            fprintf([ ...
                '\nReference simulations did not converge. ' ...
                'Skipping individual model simulations.\n' ...
                ]);
            caseCount = caseCount + length(parOut);
        else
            %% Run Loop
            for iIn = 1:length(parIn)

                switch CONV_STUDY_TYPE
                    case "space"
                        ih   = iOut;
                        iSeg = iIn;

                        if ~isnan(nSegRefError)
                            nSegRefErrorCur = nSegRefError;
                        else
                            nSegRefErrorCur = nSeg(iSeg);
                        end
                        hRefErrorCur = h(ih);
                    case "time"
                        ih   = iIn;
                        iSeg = iOut;

                        if ~isnan(hRefError)
                            hRefErrorCur = hRefError;
                        else
                            hRefErrorCur = h(ih);
                        end
                        nSegRefErrorCur = nSeg(iSeg);
                end

                fprintf('\n\nStarting Case %d/%d/%d: Mat=%d h=%.2E=2^%.2f nSeg=%d \n', ...
                    iMat, iOut, iIn, iMat, h(ih), log2(h(ih)), nSeg(iSeg));

                % Assign Input Parameters for Current Simulation Case
                simPars.h = h(ih);

                caseCount = caseCount + 1;
                fprintf('     Overall: Case %d / %d\n\n', ...
                    caseCount, length(beamPars)*length(h)*length(nSeg) )

                % Compute Reference and Initial Configuration
                [simPars.gRef, simPars.xiRef] = beamSimRefConf(nSeg(iSeg), beamPars(iMat));
                %[simPars.g0, simPars.xi0]     = beamSimInitialConf(nSeg(iSeg), beamPars(iMat));
                simPars.g0  = simPars.gRef;
                simPars.xi0 = simPars.xiRef;

                checkCFLLimit(beamPars(iMat), nSeg(iSeg), 'h', h(ih));


                %% Simulate Models

                for iSim = 1:length(beamModels)

                    if ~DEBUG_SCRIPT
                        close all
                    end

                    %% Simulate Model and Compute Error

                    fprintf('\nSimulating Model %d/%d...\n', iSim, length(beamModels))

                    % Create model object
                    beamSim = beamSimulation;
                    beamSim.simPars  = simPars;
                    beamSim.beamPars = beamPars(iMat);
                    beamSim.simModel = beamModels(iSim);

                    % Add external forcing at the tip
                    % (Smooth impulse (cosine) with specified max. force
                    % and end time)
                    beamSim.simPars.force_scaling_mode = 1;
                    beamSim.simPars.f_node_s           = [zeros(6,nSeg(iSeg)),fMax];
                    beamSim.simPars.force_tEnd         = fTEnd;

                    % Simulate model
                    beamSim = beamSim.simulateModel;

                    if beamSim.simRes.metaDataSim.exitCode
                        for iRefSim = 1:length(refBeamModels)
                            % Compute simulation errors
                            if beamSimRef(iRefSim).simRes.metaDataSim.exitCode
                                [simErrors(iMat, ih, iSeg, iSim, iRefSim), ~] = simResComputeReferenceErrors( ...
                                    hRefErrorCur, nSegRefErrorCur, beamPars(iMat).L, ...
                                    beamSim.simRes.simData, beamSimRef(iRefSim).simRes.simData ...
                                    );
                                if isnan(simErrors(iMat, ih, iSeg, iSim, iRefSim))
                                    warning("Error while computing simulation error: NaN returned.")
                                end
                            end
                        end

                        % Store simulation time
                        simTimeTotal(iMat, ih, iSeg, iSim) = beamSim.simRes.metaDataSim.TotalTime;
                    end


                    %% Plot results (if enabled)
                    if PLOT_RESULTS_INDIV
                        disp('Generating Individual Plots...')

                        % General results
                        beamSim.plotSimResults( ...
                            "savePlotsJPEG", true, ...
                            "savePlotsFig", true, ...
                            "saveDir", saveDirCase(iMat, ih, iSeg), ...
                            "computeEnergy", true, ...
                            "hPlot", max([h(ih), 2^-13])...
                            );

                        % Comparison with reference simulation
                        figHandlesRefComp = beamSim.plotReferenceComparison( ...
                            beamSimRef(1), ...
                            "savePlotsJPEG", true, ...
                            "savePlotsFig", false, ...
                            "saveDir", saveDirCase(iMat, ih, iSeg) ...
                            );

                        drawnow;
                    end

                    %% Debug Plots

                    if DEBUG_SCRIPT
                        % Velocity Comparison
                        plotColors6 = lines(6);
                        fh = figure;
                        ax = axes(fh);
                        plot(beamSim.simRes.simData.tout, squeeze(beamSim.simRes.simData.eta(:,end,:)), '-o');
                        hold on
                        plot(beamSimRef.simRes.simData.tout, squeeze(beamSimRef.simRes.simData.eta(:,end,:)), '--x');
                        legend(...
                            [cellstr(num2str((1:6).', '$\\omega_%d$ Ref.'));...
                            cellstr(num2str((1:6).', '$\\omega_%d$ Sim.'))],...
                            'Interpreter', 'latex' ...
                            );
                        colororder(ax, plotColors6);

                        % Solver Stats
                        nameStr = 'Case Sim ';
                        if ~isempty(beamSim.simRes.metaDataSteps.ImplicitIterations)
                            tout = beamSim.simRes.simData.tout;

                            figHandles.SolverStats = figure( ...
                                'Name', [nameStr, 'Solver Stats'], 'NumberTitle','off');

                            t = tiledlayout(figHandles.SolverStats, 2,2);

                            % Residual Error
                            ax = nexttile(t);
                            plot(ax, tout, beamSim.simRes.metaDataSteps.ImplicitError');
                            title(ax, 'Residual Error', 'interpreter', 'latex')
                            xlim(ax, [0, tout(end)])
                            grid on

                            % Implicit solver error flag
                            ax = nexttile(t);
                            plot(ax, tout, beamSim.simRes.metaDataSteps.ExitFlag', '-o', 'MarkerSize', 3)
                            title(ax, 'Implicit Solver Error Flag', 'interpreter', 'latex')
                            xlim(ax, [0, tout(end)])
                            grid on

                            % Nr. of iterations (of the implicit solver)
                            ax = nexttile(t, [1,2]);
                            plot(ax, tout, beamSim.simRes.metaDataSteps.ImplicitIterations', '-o');
                            title(ax, 'Nr. of Iterations', 'interpreter', 'latex')
                            xlim(ax, [0, tout(end)])
                            grid on

                            xlabel(ax, 'time $t$ / s', 'interpreter', 'latex')
                        end
                    end
                    %%


                    %% Save data (with or without full simulation data)
                    disp('Saving Individual Data...');

                    % Remove the detailed simulation results data;
                    % only keeps lightweight aggregate simulation metadata
                    if ~SAVE_DATA
                        beamSim.clearSimData();
                    end

                    % Save results
                    beamSim.saveObject( ...
                        fullfile( ...
                        saveDirCase(iMat, ih, iSeg), ...
                        strcat("beamMdl ", beamSim.simModel.modelName, ".mat") ) ...
                        );

                end % for iSim
            end % for iIn


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
                        beamSimRef(1).simRes.simData, ...
                        beamPars(iMat), true, ...
                        movieSavePath ...
                        );
                catch ME
                    warning(ME.identifier, ...
                        'Error during animation:\n %s', ME.message);
                end
            end


            %% Plot and save reference simulations
            % (Do this at the end since we clear the simulation results data
            % from the variable, which we need however to compute the
            % simulation error above
            for iRefSim = 1:length(refBeamModels)
                % Plot results (if enabled)
                if PLOT_RESULTS_INDIV
                    disp('Generating Individual Plots...')

                    % General results
                    beamSimRef(iRefSim).plotSimResults( ...
                        "savePlotsJPEG", true, ...
                        "savePlotsFig", true, ...
                        "saveDir", saveDirCase(iMat, ih, iSeg), ...
                        "computeEnergy", true, ...
                        "hPlot", 2^-13......
                        );
                end

                % Save data (with or without full simulation data)
                disp('Saving Individual Data...');

                % Save results;
                % Always clear simulation results data since it will be
                % usually ver large
                % ToDo: Alternatively, resample at lower space/time grid
                beamSimRef(iRefSim).clearSimData();
                beamSimRef(iRefSim).saveObject( ...
                    fullfile( ...
                    saveDirCase(iMat, ih, iSeg), ...
                    strcat("beamMdl ", beamSimRef(iRefSim).simModel.modelName, ".mat") ) ...
                    );
            end

        end % if refMdl Sim successful

    end % for iOut
end % for iMat

%% Save overall data

disp('Saving Overall Results Data...')

try
    save( ...
        fullfile(saveDirAll, "convStudyResults.mat"), ...
        'simErrors','simTimeTotal', ...
        'beamPars', 'nSeg', 'h', 'beamModels', 'refBeamModels',...
        '-mat', '-v7.3' ...
        );
catch ME
    warning(ME.identifier, 'Could not save data:\n %s', ME.message);
end


%% End script
simStopTime = datetime;

fprintf('Output folder: %s\n', saveDirAll);

fprintf(...
    'Finished. Time: %s, Total duration: %s (hrs/min/s)\n', ...
    string(simStopTime, 'dd.MM.yy, HH:mm:ss'), ...
    string( duration(simStopTime-simStartTime, 'Format', 'hh:mm:ss') ) ...
    );

% Turn off diary
diary('off')
