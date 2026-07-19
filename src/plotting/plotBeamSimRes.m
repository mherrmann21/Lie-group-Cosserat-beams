function figHandles = plotBeamSimRes(simRes, xiRef, modelName)
    %% Plot Simulation Results of a Beam Simulation
    %
    % Plots the results of a beam simulation.
    % The simulation results are stored in the simRes struct with the
    % standardized fields.
    %
    % Output: Struct with the figure handles of all generated plots (which
    % can be used e.g., to save the plots later).
    %
    % Maximilian Herrmann
    % Chair of Automatic Control
    % TUM School of Engineering and Design
    % Technical University of Munich

    arguments
        % Standard simRes struct
        simRes  (1,1) beamSimRes

        % 2D Array of reference deformations with dimensions (6, nSeg)
        xiRef   (6,:) double

        % Model Name
        modelName (1,1) string
    end

    %% Preparation

    simData = simRes.simData;

    % Get nr. of nodes and segments
    nNodes = size(simData.R, 3);
    nSeg   = nNodes - 1;

   

    % Time vector
    tout = simData.tout;

    % Generate string for figure names
    if ~isempty(modelName)
        nameStr = sprintf('%s: ', modelName);
    else
        nameStr = '';
    end


    %% Plot node data
    figHandles = plotBeamNodeData(simData, xiRef, "name",nameStr);



    %% Orthogonality error (all in one plot)
    if ~isempty(simRes.metaDataSteps.orthErrorR)
        figHandles.orthError = figure( ...
            'Name', [nameStr, 'OrthError'], 'NumberTitle','off');

        ax = axes(figHandles.orthError);

        plot(ax, tout, simRes.metaDataSteps.orthErrorR');
        title('Orthogonality Error $\left \| I - R^T R \right \|$', 'interpreter', 'latex')
        grid on
        xlabel('time $t$ / s', 'interpreter', 'latex')
        ylabel('$\left \| I - R^T R \right \|$', 'interpreter', 'latex')

    end


    %% Energies
    if ~isempty(simRes.E.T)
        figHandles.Energies = figure( ...
            'Name', [nameStr, 'Energies'], 'NumberTitle','off');

        ax = axes(figHandles.Energies);

        plot(ax, tout, simRes.E.H, 'LineWidth', 2);
        hold on;
        plot(ax, tout, simRes.E.T);
        plot(ax, tout, simRes.E.U);
        plot(ax, tout, simRes.E.V);

        title(ax, 'Energies', 'interpreter', 'latex')
        grid on
        xlabel(ax, 'time $t$ / s', 'interpreter', 'latex')
        ylabel(ax, 'Energy / J', 'interpreter', 'latex')

        legend(ax, 'Total $H$', 'Kinetic $T$', 'Potential $U$', 'Strain $V$', ...
            'interpreter', 'latex');
    else
        % Still create figure handle struct field so that the output struct
        % always has the same fields
        figHandles.Energies = [];
    end


    %% Solver Statistics
    if ~isempty(simRes.metaDataSteps.ImplicitIterations)

        figHandles.SolverStats = figure( ...
            'Name', [nameStr, 'Solver Stats'], 'NumberTitle','off');

        t = tiledlayout(figHandles.SolverStats, 2,1);

        % Residual Error
        ax = nexttile(t);
        plot(ax, tout, simRes.metaDataSteps.ImplicitError');
        title(ax, 'Residual Error', 'interpreter', 'latex')
        xlim(ax, [0, tout(end)])
        grid on
        xlabel(ax, 'time $t$ / s', 'interpreter', 'latex')

        % Nr. of iterations and error flag of the implicit solver
        ax = nexttile(t);
        yyaxis left
        plot(ax, tout, simRes.metaDataSteps.ImplicitIterations', '-o');
        ylabel('Nr. of iterations', 'Interpreter','latex');
        
        yyaxis right
        plot(ax, tout, simRes.metaDataSteps.ExitFlag', '-o', 'MarkerSize', 3)
        ylabel('Error Flag', 'Interpreter','latex');

        title(ax, 'Nr. of Iteration and Error Flag', 'interpreter', 'latex')
        xlim(ax, [0, tout(end)])
        grid on

        xlabel(ax, 'time $t$ / s', 'interpreter', 'latex')
    else
        % Still create figure handle struct field so that the output struct
        % always has the same fields
        figHandles.SolverStats = [];
    end



    %% 3D Beam Configuration snapshots

    % Snapshot times
    hSnapShot = tout(end) / 25;

    % Constant snapshot times
    %hSnapShot = 0.05;

    tQuery = 0:hSnapShot:simRes.simData.tout(end);

    % Interpolate configuration in time
    gInput = SE3Matrix(simRes.simData.R,simRes.simData.x);
    gQuery = interpn( ...
        1:4, 1:4, 1:size(gInput,3), simRes.simData.tout, ...
        gInput, ...
        1:4, 1:4, 1:size(gInput,3), tQuery');

    snapShotColormap = winter(size(gQuery,4));

    %%% Actual Plot
    figHandles.beamSnapshots = figure( ...
        'Name', [nameStr, 'snapshots'], 'NumberTitle','off');
    axSnapShots = axes(figHandles.beamSnapshots);
    hold on;
    title( ...
        axSnapShots, sprintf('Configuration Snapshots (%.3f s steps)', hSnapShot), ...
        'Interpreter','latex');

    %%% Plot individual snaphshot configurations
    % Disable warnings that may appear for failed simulations
    warning('off', 'MATLAB:illConditionedMatrix');
    warning('off', 'MATLAB:nearlySingularMatrix');
    warning('off', 'MATLAB:singularMatrix');

    for iStep = 1:size(gQuery,4)
        elasticBeamSimple(gQuery(:,:,:,iStep), ...
            'name', sprintf('$t = %.2f$ s', tQuery(iStep)), ...
            'color', snapShotColormap(iStep,:), 'axis', axSnapShots);
    end

    %legend('Interpreter','latex');
    axis equal
    grid on

    width = 0.015;
    xlim([ ...
        min(gQuery(1,4,:,:),[],'all')-width, ...
        max(gQuery(1,4,:,:),[],'all')+width ...
        ]);
    ylim([ ...
        min(gQuery(2,4,:,:),[],'all')-width, ...
        max(gQuery(2,4,:,:),[],'all')+width ...
        ]);
    zlim([ ...
        min(gQuery(3,4,:,:),[],'all')-width, ...
        max(gQuery(3,4,:,:),[],'all')+width ...
        ]);

    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    zlabel('$z$','Interpreter','latex');

    view([37.5, 30]);

    axSnapShots.Colormap = winter;
    ch = colorbar;
    ch.Label.Interpreter = 'latex';
    ch.TickLabelInterpreter = 'latex';
    %ch.FontSize = 12;
    clim(axSnapShots, [tout(1), tout(end)]);
    ylabel(ch, 'time in s', 'Interpreter', 'latex','FontSize',12);



    %% Surf Plot Xi
    figHandles.surfStrains = figure( ...
        'Name', [nameStr, 'Surf Strains'], 'NumberTitle','off');

    tlXi = tiledlayout(figHandles.surfStrains, 'flow');
    xiRefMat = repmat(xiRef, [1,1, length(tout)]);
    for iVal = 1:6
        nexttile(tlXi);
        axis equal
        surf( ...
            tout, linspace(0,1,nSeg), ...
            squeeze(simData.xi(iVal,:,:)-xiRefMat(iVal,:,:)), ...
            'edgecolor','none', ...
            'FaceColor', 'interp'...
            );
        title(sprintf('Deformation $(\\xi_%d-\\bar{\\xi}_%d)$', iVal,iVal), 'Interpreter', 'latex');
        colormap jet

        xlim([tout(1), tout(end)]);
        xlabel('time $t$ / s', 'interpreter', 'latex');
        ylim([0, 1]);
        ylabel('norm. beam length', 'interpreter', 'latex');
        zlabel(sprintf('$(\\xi_%d-\\bar{\\xi}_%d)$', iVal,iVal), 'interpreter', 'latex');
        view([-135, 25]);
    end


    %% Surf Plot eta
    figHandles.surfVel = figure( ...
        'Name', [nameStr, 'Surf Velocities'], 'NumberTitle','off');
    tlEta= tiledlayout(figHandles.surfVel, 'flow');
    for iVal = 1:6
        nexttile(tlEta);
        axis equal
        surf( ...
            tout, linspace(0,1,nNodes), ...
            squeeze(simData.eta(iVal,:,:)), ...
            'edgecolor', 'none', ...
            'FaceColor', 'interp'...
            );
        title(sprintf('Velocity $\\eta_%d$', iVal), 'Interpreter', 'latex');
        colormap jet

        xlim([tout(1), tout(end)]);
        xlabel('time $t$ / s', 'interpreter', 'latex');
        ylim([0, 1]);
        ylabel('norm. beam length', 'interpreter', 'latex');
        zlabel(sprintf('$\\eta_%d$', iVal), 'interpreter', 'latex');
        view([-135, 25]);
    end

    %% Surf Plot x
    figHandles.surfPos = figure( ...
        'Name', [nameStr, 'Surf Positions'], 'NumberTitle','off');
    tlX = tiledlayout(figHandles.surfPos, "flow");
    varStrings = ["x","y","z"];
    for iVal = 1:3
        nexttile(tlX);
        axis equal
        surf( ...
            tout, linspace(0,1,nNodes), ...
            squeeze(simData.x(iVal,:,:)), ...
            'edgecolor','none', ...
            'FaceColor', 'interp'...
            );
        title(strcat("Position $", varStrings(iVal), "$"), 'Interpreter', 'latex');
        colormap jet

        xlim([tout(1), tout(end)]);
        xlabel('time $t$ / s', 'interpreter', 'latex');
        ylim([0, 1]);
        ylabel('norm. beam length', 'interpreter', 'latex');
        zlabel(strcat("$", varStrings(iVal), "$ / m"), 'interpreter', 'latex');
        view([-135, 25]);
    end




end
