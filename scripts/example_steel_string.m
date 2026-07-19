%% Example Simulation: Dynamic Simulation of a Steel String 
% Used integrator: Relative-kinematic LGVI with Broyden solver
% Used beam model: Inextensible Kirchhoff beam (bending and torsion only)
%
% Note that the steel string results in a fairly stiff model due to its
% slenderness and material properties; hence, it is advantageous to
% simulate it as a reduced Kirchhoff beam (without the stiff torsion and
% extension modes), which results in much faster computation times.
% Simulations with a Simo-Reissner beam would require much finer time
% steps.
%
% Maximilian Herrmann
% Chair of Automatic Control
% TUM School of Engineering and Design
% Technical University of Munich

clear
close all

%% Script Settings

% Save plots and animations (if generated) and beam simulation object?
SAVE_RESULTS = false;

% Save simulation data?
SAVE_DATA = false;

% Generate animation (and save, if SAVE_RESULTS is true)?
ANIMATE_RESULTS = 1;

% Directory where the results subfolder will be created
RESULTS_DIR = fullfile(getRootFolder, "results", "simulations");


%% Model Parameters

% Nr. of segments n
nSeg = 2^3;

%% Beam Parameters

%beamPars = beamParams_LLA11_rubberRod();
%beamPars = beamParams_slenderRubberRod();
beamPars = beamParams_LLA11_steelString;
%beamPars = beamParams_mbsd_soft_rod;

% Explicitly set damping, if desired
% From LLA12
beamPars.d = [
    2e-4
    2e-4
    8e-6
    1e-1
    1e-1
    2e-0
    ];


%% Compute Beam Reference Configuration
% Configuration of the undeformed beam

simPars = beamSimPars;
[simPars.gRef, simPars.xiRef] = beamSimRefConf(nSeg, beamPars);

% Visualize reference configuration
visualizeBeamConfig(simPars.gRef, beamPars, 'Reference Configuration');


%% Compute Initial conditions

% Start simulation in reference configuration
simPars.g0  = simPars.gRef;
simPars.xi0 = simPars.xiRef;

% Visualize initial configuration
visualizeBeamConfig(simPars.g0, beamPars, 'Initial Configuration');


%% Simulation settings

% Simulation parameters
simPars.tEnd = 6;

% Create simulation object
beamSim = beamSimulation;
beamSim.simPars  = simPars;
beamSim.beamPars = beamPars;


%% Simulate as Kirchhoff Beam with Relative-Kinematic Broyden LGVI

% Define beam model to simulate
beamSim.simModel.modelName        = "RelKin KH Broyden";
beamSim.simModel.reducedParams.Ba = [ eye(3); zeros(3)];
beamSim.simModel.reducedParams.Bc = [ zeros(3); eye(3)];
beamSim.simModel.funHandle        = @(beamPars, simPars, simModel) ...
    beamMdlRelKinVarInt_Broyden_mex(simPars, beamPars, simModel.solverConfig, simModel.reducedParams.Ba, simModel.reducedParams.Bc);

% Solver settings
beamSim.simModel.solverConfig.errorMargin = 1e-10;
beamSim.simModel.solverConfig.maxIterations = 25;
beamSim.simModel.solverConfig.UseExactJabocobian = true;
beamSim.simModel.solverConfig.JacobianIterationThreshold = 5;

% Set time step
beamSim.simPars.h = 2^-8.7;
checkCFLLimit(beamPars, nSeg, 'h', beamSim.simPars.h);

% Simulate model
beamSim = beamSim.simulateModel("consoleOutput", true, "accurateTiming", false);


%% Raw Solver Stats (non-interpolated)

figHandles.SolverStats = figure( ...
    'Name', 'Solver Stats', 'NumberTitle','off');

t = tiledlayout(figHandles.SolverStats, 2,1);

% Implicit solver error flag
ax = nexttile(t);
plot(ax, beamSim.simRes.metaDataSteps.ExitFlag', '-o', 'MarkerSize', 3)
title(ax, 'Implicit Solver Error Flag', 'interpreter', 'latex')
xlim(ax, [0, length(beamSim.simRes.simData.tout)])
grid on

% Residual Error
ax = nexttile(t);
yyaxis left
plot(ax, beamSim.simRes.metaDataSteps.ImplicitError');
title(ax, 'Residual Error', 'interpreter', 'latex')
xlim(ax, [0, length(beamSim.simRes.simData.tout)])
grid on
ylabel('Residual', 'Interpreter','latex');

% Nr. of iterations (of the implicit solver)
%ax = nexttile(t, [1,2]);
yyaxis right
plot(ax, beamSim.simRes.metaDataSteps.ImplicitIterations', '-o');
title(ax, 'Residual and Nr. of Iterations', 'interpreter', 'latex')
xlim(ax, [0, length(beamSim.simRes.simData.tout)])
grid on
ylabel('Nr. of Iterations', 'Interpreter','latex');

xlabel('Simulation Steps', 'interpreter', 'latex')


%% Plot results 

% Create output folder
subFolder = sprintf("%s_steel_string_n=%d_h=2^%.2f", ...
    string(datetime, 'yyMMdd_HHmm'), nSeg, log2(beamSim.simPars.h) ...
    );
saveDir = fullfile(RESULTS_DIR, subFolder);
if (SAVE_RESULTS || SAVE_DATA) && ~isfolder( saveDir )
    mkdir( saveDir ); %#ok<*UNRCH>
end

% Only plot if simulation was successful
    if beamSim.simRes.metaDataSim.exitCode

        disp('Generating Individual Plots...')

        beamSim.plotSimResults( ...
            "savePlotsJPEG", SAVE_RESULTS, ...
            "savePlotsFig", SAVE_RESULTS, ...
            "saveDir", saveDir, ...
            "computeEnergy", true, "hPlot", 2^-12 ...
            );
        drawnow;
    end



%% Animate Beam Results for Visual Check

if ANIMATE_RESULTS && beamSim.simRes.metaDataSim.exitCode
    disp('Animating Results...')

    movieSavePath = fullfile(saveDir,"animation");

    try
        animateBeamSimRes( ...
            beamSim.simRes.simData, ...
            beamPars, SAVE_RESULTS, ...
            movieSavePath ...
            );
    catch ME
        warning(ME.identifier, ...
            'Error during animation:\n %s', ME.message);
    end
end


%% Save data (with or without full simulation data)

if SAVE_RESULTS
    disp('Saving Individual Data...');

    % Remove the detailed simulation results data;
    % only keeps "lightweight" overall simulation meta data
    if ~SAVE_DATA
        beamSim.clearSimData();
    end

    % Save results
    beamSim.saveObject( ...
        fullfile( ...
        saveDir, ...
        strcat("beamMdl ", beamSim.simModel.modelName, ".mat") ) ...
        );

end


%% End script

disp('Finished.')
