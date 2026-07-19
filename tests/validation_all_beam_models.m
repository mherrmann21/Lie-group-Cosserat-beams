%% Validate all beam models: Run simulation for all models
% Simulate all models with the same simulation conditions and plot the
% results for comparison / validation.
%
% Maximilian Herrmann
% Chair of Automatic Control
% TUM School of Engineering and Design
% Technical University of Munich


%% Preparation

clear
close all


%% Script Settings

PLOT_RESULTS_INDIV = 1;


%% Model Parameters

% Nr. of segments n
nSeg = 2^4;


%% Beam Parameters

%beamPars = beamParams_default();
beamPars = beamParams_slenderSteelRod('E', 1e9);

% Explicitly set damping, if desired
beamPars.d = ones(6,1)*1e-4;


%% Reference and Initial Configuration

% Compute Beam Reference Configuration
[gRef, xiRef] = beamSimRefConf(nSeg, beamPars);

% Visualize reference configuration
visualizeBeamConfig(gRef, beamPars, 'Reference Configuration');

% Compute Initial conditions

[g0, xi0] = beamSimInitialConf(nSeg, beamPars);

% Visualize initial configuration
visualizeBeamConfig(g0, beamPars, 'Initial Configuration');


%% Get beam models to compare / validate

% Solver settings
solverConfig = beamSolverConfig;
solverConfig.errorMargin   = 1e-8;
solverConfig.maxIterations = 200;
solverConfig.JacobianIterationThreshold = 5;

% Simulation parameters
simPars = beamSimPars;
simPars.g0    = g0;
simPars.gRef  = gRef;
simPars.xiRef = xiRef;
simPars.xi0   = xi0;

simPars.tEnd  = 0.25;%1;
simPars.h     = 2^-16;

simPars.f_node_s = zeros(6,nSeg+1);
simPars.f_node_s(:,end) = [0 0 0 0 0 1]*1e-1;

beamModels = defineSimStudyBeamModels('solverConfig', solverConfig);

% Add continuous-time model
beamModels(end+1).modelName      = "AbsKin ode15s SR";
beamModels(end).solverConfig     = solverConfig;
beamModels(end).reducedParams.Ba = [ eye(3); zeros(3)];
beamModels(end).reducedParams.Bc = [ zeros(3); eye(3)];
beamModels(end).funHandle        = @(beamPars, simPars, simModel) ...
    simulateBeam_absKin_cont(beamPars, simPars, @ode15s);


%% Simulate Models

for iSim = 1:length(beamModels)

    fprintf('\nSimulating Model %d/%d...\n', iSim, length(beamModels))

    % Create model object
    beamSim = beamSimulation;
    beamSim.simPars  = simPars;
    beamSim.beamPars = beamPars;
    beamSim.simModel = beamModels(iSim);

    % Simulate model
    beamSim = beamSim.simulateModel;

    % Compute energy behavior (to provide qualitative
    % insight into simulation accuracy)
    %beamMdl.computeEnergyEvolution();


    %% Plot results (if enabled)
    if PLOT_RESULTS_INDIV
        disp('Generating Individual Plots...')

        % General results
        beamSim.plotSimResults("hPlot", max([2^-9, simPars.h]));

        drawnow;
    end

end

disp('Finished.')
