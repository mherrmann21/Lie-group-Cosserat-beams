%% Example: Dynamic Simulation of a Rubber Rod using the absolute-kinematic ODE model
% For illustration, the simulation is done with several beam models:
%  1. Relative-kinematic Kirchhoff beam model with LGVI and Broyden solver
%  2. Absolute-kinematic Simo-Reissner beam model with LGVI
%  3. Absolute-kinematic Simo-Reissner beam model with ODE solver (ode15s)
%
% The simulated beam is fairly soft, so the numerical models have
% comparatively low stiffness and can be simulated efficiently even with the
% Simo-Reissner models.
%
% Maximilian Herrmann
% Chair of Automatic Control
% TUM School of Engineering and Design
% Technical University of Munich

clear
close all

%% Script Settings

% Generate animation?
ANIMATE_RESULTS = true;

% Save animation?
SAVE_MOVIE = true;


%% Model Parameters

% Nr. of segments n
nSeg = 10;


%% Beam Parameters

beamPars = beamParams_LLA11_rubberRod;

% Explicitly set damping, if desired
beamPars.d = ones(6,1)*1e-4;


%% Compute Beam Reference Configuration
% Configuration of the undeformed beam

[gRef, xiRef] = beamSimRefConf(nSeg, beamPars);

% Visualize reference configuration
visualizeBeamConfig(gRef, beamPars, 'Reference Configuration');


%% Compute Initial conditions

[g0, xi0] = beamSimInitialConf(nSeg, beamPars);

% Visualize initial configuration
visualizeBeamConfig(g0, beamPars, 'Initial Configuration');


%% Simulation parameters
simPars = beamSimPars;
simPars.gRef = gRef;
simPars.g0 = g0;
simPars.xiRef = xiRef;
simPars.tEnd = 1;

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
beamSim.simModel.solverConfig.JacobianIterationThreshold = 4;

% Set time step
beamSim.simPars.h = 2^-10;
checkCFLLimit(beamPars, nSeg, 'h', beamSim.simPars.h);


% Simulate model
beamSim = beamSim.simulateModel("consoleOutput", true, "accurateTiming", false);


%% Simulate as Simo-Reissner Beam with Absolute-Kinematic LGVI

% Define beam model to simulate
beamSim.simModel.modelName        = "AbsKin SR General";
beamSim.simModel.funHandle        = @(beamPars, simPars, simModel) ...
    beamMdlAbsKinLGVI_general_mex(simPars, beamPars, simModel.solverConfig, zeros(6,0), zeros(6,0));

% Solver settings
beamSim.simModel.solverConfig.errorMargin = 1e-8;
beamSim.simModel.solverConfig.maxIterations = 25;
beamSim.simModel.solverConfig.UseExactJacobian = true;
beamSim.simModel.solverConfig.JacobianIterationThreshold = 4;

% Set time step (must be small for SR model)
beamSim.simPars.h = 2^-14;
checkCFLLimit(beamPars, nSeg, 'h', beamSim.simPars.h);

% Simulate model
beamSim = beamSim.simulateModel("consoleOutput", true, "accurateTiming", false);


%% Integration with ODE solver

disp('Starting integration...')

% ODE integrator options
opts = odeset();

% Define beam model/integrator
beamSim.simModel.modelName        = "AbsKin ode15s SR";
beamSim.simModel.reducedParams.Ba = [ eye(3); zeros(3)];
beamSim.simModel.reducedParams.Bc = [ zeros(3); eye(3)];
beamSim.simModel.funHandle        = @(beamPars, simPars, simModel) ...
    simulateBeam_absKin_cont(beamPars, simPars, @ode15s, opts);

% Simulate model
beamSim = beamSim.simulateModel("consoleOutput", true, "accurateTiming", false);


%% Plot Results

disp('Generating Plots...')
beamSim.plotSimResults;


%% Animate Results
% And save to video

if ANIMATE_RESULTS
    disp('Animating Results...')
    movieFileName = sprintf("%s_simulation_rubber_rod", string(datetime, 'yyMMdd_HHmm'));
    movieSavePath = fullfile(getRootFolder, "results", "simulations", movieFileName);
    animateBeamSimRes(beamSim.simRes.simData, beamPars, SAVE_MOVIE, movieSavePath);
end

disp('Finished.')
