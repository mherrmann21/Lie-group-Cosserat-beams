%% Validate the external forces time step function
%
% Maximilian Herrmann
% Chair of Automatic Control
% TUM School of Engineering and Design
% Technical University of Munich

clear
close all


%% Test parameters
nNodes = 9;
h = 2^-6;
tEnd = 10;

%% Define Force

% Time at which the force interval ends
tForceEnd = 5;

%%% Example for force input as time-value sample pairs
% Nr. of samples in the transient part
nSampleRise = round(tForceEnd/2e-6);

scaleVec = [0,(1-cos(pi*linspace(0,1,nSampleRise)))/2, 1];
tVec = [linspace(0,tForceEnd,nSampleRise+1), 3000];

fConst = [0 0 0 0 0 600 ]';

simPars = beamSimPars;
simPars.f_node_s = [zeros(6,nNodes-1),fConst];
simPars.f_node_tVec = tVec;
simPars.f_node_sVec = scaleVec;
simPars.force_tEnd = tForceEnd;


%% Compute test

tout = 0:h:tEnd;
simPars.g0 = zeros(4,4,nNodes);
fCalc_s = zeros(6,nNodes,length(tout));

for iMode = 0:4

    fprintf('Force Mode = %d\n', iMode);
    simPars.force_scaling_mode = iMode;

    tic
    for iStep = 1:length(tout)
        [f_node_k_b, fCalc_s(:,:,iStep)] = getExtStepNodeForces(simPars, tout(iStep));
    end
    toc

    %% Plot forces at beam end
    figure('Name', sprintf('Mode %d Ext. Forces Tip', iMode), 'NumberTitle','off');
    plot(tout,  squeeze(fCalc_s(:,end,:)) );
    title('External Tip Forces over Time', 'Interpreter','latex');
    xlabel('time / s', 'Interpreter','latex');
    ylabel('force components', 'Interpreter','latex');
    grid on;
    legend(string(1:6));
    drawnow;


    %% Surf Plot

    fh = figure( 'Name', sprintf('Mode %d Ext. Forces Surf', iMode), 'NumberTitle','off');
    tlEta= tiledlayout(fh, 'flow');
    for iVal = 1:6
        nexttile(tlEta);
        axis equal
        surf( ...
            tout, linspace(0,1,nNodes), ...
            squeeze(fCalc_s(iVal,:,:)), ...
            'edgecolor', 'none', ...
            'FaceColor', 'interp'...
            );
        title(sprintf('Ext. Force $f_%d$', iVal), 'Interpreter', 'latex');
        colormap jet

        xlim([tout(1), tout(end)]);
        xlabel('time $t$ / s', 'interpreter', 'latex');
        ylim([0, 1]);
        ylabel('norm. beam length', 'interpreter', 'latex');
        zlabel(sprintf('$f_%d$', iVal), 'interpreter', 'latex');
        view([-135, 25]);
    end
    drawnow;

end