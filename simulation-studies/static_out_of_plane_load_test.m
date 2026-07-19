%% Static Test for RelKin Paper: Out-of-Plane Load
%
% Maximilian Herrmann
% Chair of Automatic Control
% TUM School of Engineering and Design
% Technical University of Munich

clear
close all

%% Script Settings

SOL_STAT = 1;
SOL_DYN  = 1;

SAVE_RESULTS = 1;

RESULTS_DIR = fullfile(getRootFolder, "results", "runs");

% Use accurate timing (with timeit) for dynamic simulations?
ACCURATE_TIMING = false;

% Use exponential map instead of Cayley map?
USE_EXP = false;


%% Simulation parameters

% Segment numbers to investigate
nSegVec = [8, 32];

% Beam parameters
beamPars = beamParams_outOfPlaneLoad();


%% Prepare Output Folder and Start Logging

subFolderName = sprintf( "%s_mbsd_test_outOfPlaneLoad_useExp=%s", string(datetime, 'yyMMdd_HHmm'), string(USE_EXP));
saveDir = fullfile(RESULTS_DIR, subFolderName);

simStartTime = datetime;
if SAVE_RESULTS
    mkdir( saveDir );

    diary(fullfile(saveDir, 'simulation.log'));
    fprintf('Starting Log. Time: %s\n', string(simStartTime, 'dd.MM.yy, HH:mm:ss'));
end

%% Simulate cases

nSims = 5; % 2 static + 3 dynamic simulations

endPos = nan(3, nSims, length(nSegVec));

for iSeg = 1:length(nSegVec)
    nSeg = nSegVec(iSeg);
    fprintf("\n Starting case nSeg = %d (case %d/%d)...\n\n", nSeg, iSeg, length(nSegVec));

    %% Initial configuration
    % Beam is aligned with the global X axis

    % Segment length
    l = beamPars.L / nSeg;

    % Configuration of first node
    RRef0 = eul2rotm([pi/2, pi/2, 0]);
    g0 = SE3Matrix(RRef0, zeros(3,1));

    % Vector of constant deformations used for all segments
    xiConst = [2*pi/beamPars.L/8; 0; 0; 0; 0; 1];

    % Array of segment deformations
    xiRef = repmat(xiConst, [1, nSeg]);

    % Forward kinematics
    if USE_EXP
        [~, xRef, gRef] = beamRelFwdKinExp(xiRef*l, g0);
    else
        [~, xRef, gRef] = beamRelFwdKin(xiRef*l, g0);
    end

    % Visualize reference configuration
    visualizeBeamConfig(gRef, beamPars, 'Reference Configuration', ...
        'showFrames',true, 'showLabels',true);

    % Visualize the circle to make sure everything is correct
    fh2 = figure;
    ax2 = axes(fh2);
    axis equal
    hold on;
    grid on
    elasticBeamSimple(gRef, 'name', 'Ref Config', 'axis', ax2);

    th = 0:pi/50:2*pi;
    radius = 100;
    xunit = radius * cos(th) + radius;
    yunit = radius * sin(th) + 0;
    h = plot(xunit, yunit, ':');
    plot(0:100, (0:100)*(-1)+radius)


    disp('End position initial configuration:')
    disp(xRef(:,end)');


    %% Compute solution from static models

    if SOL_STAT
        simParsStat = beamSimPars;
        simParsStat.g0 = gRef;
        simParsStat.gRef = gRef;
        simParsStat.xiRef = xiRef;
        simParsStat.xi0   = xiRef;
        simParsStat.f_node_s = [zeros(6,nSeg), [0 0 0 0 0 600].'];
        simParsStat.g = 0;


        %%% Simo-Reissner
        Ba_SR = eye(6);
        Bc_SR = zeros(6,0);

        nLoadSteps = 200;
        tic;
        [gStat_SR, simData_gStat_SR, solInfo] = computeStaticEquilibriumRelKin_mex(simParsStat, beamPars, Ba_SR, Bc_SR, nLoadSteps, USE_EXP);
        tSim = toc;
        fprintf('Solved static equilibrium in %f s. Mean Iterations: %.1f,  Mean Residual: %e.\n', tSim, mean(solInfo.iterations), mean(solInfo.residual(:)));

        visualizeBeamConfig(gStat_SR, beamPars, 'Solution Static Sim. SR', ...
            'showFrames',false, 'showLabels',false);

        %elasticBeamSimple(gEqu2, 'name', sprintf('Static Simulation, nSeg=%d', nSeg2), 'color', colors(2,:), 'axis', ax);

        endPos(:,1,iSeg) = simData_gStat_SR.x(:,end);
        disp('End position final configuration:')
        disp(simData_gStat_SR.x(:,end)');


        %%% Kirchhoff
        Ba_KH = [ eye(3); zeros(3)];
        Bc_KH = [ zeros(3); eye(3)];

        nLoadSteps = 200;
        tic;
        [gStat_KH, simData_gStat_KH, solInfo] = computeStaticEquilibriumRelKin_mex(simParsStat, beamPars, Ba_KH, Bc_KH, nLoadSteps, USE_EXP);
        tSim = toc;
        fprintf('Solved static equilibrium in %f s. Mean Iterations: %.1f,  Mean Residual: %e.\n', tSim, mean(solInfo.iterations), mean(solInfo.residual(:)));

        visualizeBeamConfig(gStat_KH, beamPars, 'Solution Static Sim. KH', ...
            'showFrames',false, 'showLabels',false);

        %elasticBeamSimple(gEqu2, 'name', sprintf('Static Simulation, nSeg=%d', nSeg2), 'color', colors(2,:), 'axis', ax);

        endPos(:,2,iSeg) = simData_gStat_KH.x(:,end);
        disp('End position final configuration:')
        disp(simData_gStat_KH.x(:,end)');
    end

    %% Compute solutions from dynamic models

    if SOL_DYN
        simPars = beamSimPars;

        %% External Forces

        % Constant force at beam tip
        fConst = [0 0 0 0 0 600 ].';
        simPars.f_node_s = [zeros(6,nSeg),fConst];

        % Time at which the rise interval ends
        simPars.force_tEnd = 5;

        % Smooth rise from 0 to max.
        simPars.force_scaling_mode = 2;


        %% Simulation
        simPars.g0    = gRef;
        simPars.gRef  = gRef;
        simPars.xiRef = xiRef;
        simPars.xi0   = xiRef;
        simPars.g = 0;

        % Solver settings for dynamic simulation
        solverConfig = beamSolverConfig;
        solverConfig.errorMargin   = 5e-8;
        solverConfig.errorMarginLimit   = 1e-6;
        solverConfig.maxIterations = 200;
        solverConfig.JacobianIterationThreshold = 5;
        solverConfig.UseExactJacobian = true;

        % Get beam models to compare / validate
        beamModels = defineSimStudyBeamModels('solverConfig', solverConfig);

        % Simulation parameters
        beamPars.d = ones(6,1)*0.5e7;

        simPars.tEnd = 200;
        simPars.h    = 2^-6;

        checkCFLLimit(beamPars, nSeg, 'h', simPars.h);

        for iSim = 1:3 % Only AbsKin + RelKin Broyden (SR + KH) Models

            fprintf('\nSimulating Model %d/%d...\n', iSim, length(beamModels))

            % Create model object
            beamSim = beamSimulation();
            beamSim.simPars  = simPars;
            beamSim.beamPars = beamPars;
            beamSim.simModel = beamModels(iSim);
            beamSim.simModel.solverConfig.UseExactJacobian = true;

            %beamSim.simModel.funHandle = @(beamPars, simPars, simModel) beamMdlAbsKinLGVI_general(simPars, beamPars, simModel.solverConfig, simModel.reducedParams.Ba, simModel.reducedParams.Bc);


            % Simulate model
            beamSim = beamSim.simulateModel("accurateTiming", ACCURATE_TIMING);

            endPos(:,iSim+2,iSeg) = beamSim.simRes.simData.x(:,end,end);


            % Visualize equilibrium configuration
            gEqu_dyn = SE3Matrix( ...
                beamSim.simRes.simData.R(:,:,:,end), ...
                beamSim.simRes.simData.x(:,:,end) ...
                );

            if beamSim.simRes.metaDataSim.exitCode

                disp('Final velocities:')
                disp(beamSim.simRes.simData.eta(:,end,end-1)');

                disp('End position final configuration:')
                disp(gEqu_dyn(1:3,4,end)');

                visualizeBeamConfig(gEqu_dyn, beamPars, 'Solution Dyn. Sim');

                disp('Generating Individual Plots...')

                %beamSim.plotSimResults("hPlot", 1e-2);

                %elasticBeamSimple(gEqu_dyn2, 'name', 'Dynamic Simulation', 'color', colors(3,:), 'axis', ax);
                drawnow;

            end
        end
    end
end

%% Show results

disp(endPos)

endPosStrings = strings(nSims, length(nSegVec));
for iSeg = 1:length(nSegVec)
    for iSim = 1:nSims
        endPosStrings(iSim, iSeg) = sprintf("%.4f, %.4f, %.4f", endPos(:,iSim,iSeg));
    end
end

T = splitvars(table(endPosStrings, 'RowNames', ["static SR", "static KH", beamModels(1:3).modelName]));
T.Properties.VariableNames = string(nSegVec);
disp(T);

if SAVE_RESULTS
    writetable(T, fullfile(saveDir, "table_end_pos"),'WriteRowNames',true );
end


%% 3D Visualization
colors = tumColors();

fhVis = figure("Name", "visualization_outOfPlaneLoad");
axVis = axes(fhVis);
hold on

% Initial configuration
elasticBeam( gRef, 'width', 1, 'height', 1, ...
    "color", colors.TUMBlue1, 'showFrames', false ...
    );

% Final configuration
elasticBeam( gEqu_dyn, 'width', 1, 'height', 1, ...
    "color", colors.TUMBlue, 'showFrames', false ...
    );

% Force arrow initial conf.
quiver3( ...
    gRef(1,4,end), gRef(2,4,end), gRef(3,4,end), ...
    0, 0, 15, 'MaxHeadSize',  0.5, 'LineWidth', 1.25, 'Color', "black", ...
    'AutoScale','off' ...
    );
% Force arrow final conf.
quiver3( ...
    gEqu_dyn(1,4,end), gEqu_dyn(2,4,end), gEqu_dyn(3,4,end), ...
    0, 0, 15, 'MaxHeadSize',  0.5, 'LineWidth', 1.25, 'Color', "black", ...
    'AutoScale','off' ...
    );

% Force text
text(gRef(1,4,end)+1, gRef(2,4,end)+3, gRef(3,4,end)+15, '600~N', ...
    'Interpreter','latex', 'FontSize', 8);
text(gEqu_dyn(1,4,end)+1, gEqu_dyn(2,4,end)+3, gEqu_dyn(3,4,end)+15, '600~N', ...
    'Interpreter','latex', 'FontSize', 8);

axis equal
grid on
xlim([0,40]);
ylim([0, 80]);
zlim([-0.5, 70]);
%view([135, 25]);
view([125, 32]);

xlabel('$x$ in m', 'Interpreter', 'latex');
ylabel('$y$ in m', 'Interpreter', 'latex');
zlabel('$z$ in m', 'Interpreter', 'latex');
axVis.TickLabelInterpreter = "latex";

%%% Format and save

pdfWidth = 250; % Page width in pt
pdfAspectRatio = 1;4/3;

pdfSize = [pdfWidth, pdfWidth/pdfAspectRatio]/ 28.346; % figure size in cm

% Un-dock figure if docked
fhVis.WindowStyle = "normal";

fhVis.CurrentAxes.FontSize = 8;

% Position/Size

fhVis.Units = 'centimeters';
fhVis.Position = [20,20,pdfSize];

fhVis.PaperUnits = 'centimeters';
fhVis.PaperPosition = [0,0,pdfSize];
fhVis.PaperSize = pdfSize;

axisMarginsLeftBottom = [1.0, 0.9];
axisMarginsRightTop   = [0.3, 0.3];
axVis.Units = "centimeters";
axVis.Position = [axisMarginsLeftBottom, pdfSize-axisMarginsLeftBottom-axisMarginsRightTop];

% Save figure
if SAVE_RESULTS
    savefig(fhVis, fullfile(saveDir, fhVis.Name));
    %exportgraphics(fh, strcat(fullfile(PLOT_FOLDER, fh.Name), ".pdf"));
    print(fhVis, fullfile(saveDir, fhVis.Name), "-dpdf", '-vector');
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



%% Local functions

function [R, x, g, g_ij] = beamRelFwdKinExp(xi, g0)
    %% Compute Beam Forward Kinematics from Relative Deformations
    % i.e., compute absolute positions and rotations for each node
    %
    % Inputs:
    %  xi   Array of relative deformations with dimensions (6, nSeg)
    %
    %       Important: Must correspond to discrete *updates*, not discrete
    %       *gradients*! I.e., tau(xi) = g_a^{-1}*g_{a+1} must hold!
    %
    %  g0   SE3 configuration element of the first node (4x4 matrix)
    %
    % Outputs:
    %   R  Array of node rotation matrices with dimensions (3,3,nNodes)
    %   x  Array of node position vectors with dimensions (3,nNodes)
    %   g  Array of SE3 matrices with dimensions (4,4,nNodes) that describe
    %   the *absolute* configurations of the nodes w.r.t. g0
    %
    %   g_ij Array of SE3 matrices with dimensions (4,4,nSeg) that
    %   describe the *relative* configurations between the nodes, i.e.,
    %        g_j = g_i * g_ij
    %   holds; the definition is
    %        g_ij = tau( xi_ij)
    %   (where xi is the relative *update* as described above).
    %   Index i corresponds to the update described by xi with index i.


    % Nr. of nodes and segments
    nSeg   = size(xi, 2);
    nNodes = nSeg + 1;

    % Compute absolute configuration from relative deformations
    g    = zeros(4,4,nNodes);
    g_ij = zeros(4,4,nSeg);

    % Absolute configuration of the first node is given with g0
    g(:,:,1) = g0;

    % Compute relative update between first two nodes
    g_ij(:,:,1) = expSE3( xi(:, 1) );

    for iN = 2:nNodes
        if iN < nNodes
            g_ij(:,:,iN) = expSE3( xi(:, iN) );
        end
        g(:,:,iN) = g(:,:,iN-1) * g_ij(:,:,iN-1);
    end

    % Get arrays of rotation matrices and position vectors
    [R, x] = RxFromSE3Matrix(g);
end
