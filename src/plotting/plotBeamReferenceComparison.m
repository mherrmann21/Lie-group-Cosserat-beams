function figHandles = plotBeamReferenceComparison(simParsSim, simDataSim, simParsRef, simDataRef)
    %% Plot Comparison of a Beam Simulation with a Reference Simulation
    %
    % Inputs:
    %    Objects 'simPars' and 'simRes' of the comparison and the reference
    %    simulation.
    %
    % Outputs:
    %    Struct with figure handles to the generated figures.
    %
    % Maximilian Herrmann
    % Chair of Automatic Control
    % TUM School of Engineering and Design
    % Technical University of Munich

    arguments
        simParsSim (1,1) beamSimPars
        simDataSim (1,1) beamSimData
        
        simParsRef (1,1) beamSimPars
        simDataRef (1,1) beamSimData
    end



    %% Prepare Comparison Data
    nSeg_Sim = size(simDataSim.R, 3) - 1;
    nSeg_Ref = size(simDataRef.R, 3) - 1;

    % Segment size for unit length
    l_unit_Ref = 1/nSeg_Ref;
    l_unit_Sim = 1/nSeg_Sim;


    % Check if the sampling time and discretization step of the reference
    % simulation is a multiple of the comparison simulation
    if mod( simParsSim.h, simParsRef.h )
        error( ...
            'Sampling time of the reference simulation (%.3E) is not an integer multiple of the comparison sampling time (%.3E).', ...
            simParsRef.h, simParsSim.h ...
            );
    end
    if mod( l_unit_Sim, l_unit_Ref )
        error( ...
            'Segment length of the reference simulation (%.3E) is not an integer multiple of the comparison segment length (%.3E).', ...
            l_unit_Ref, l_unit_Sim ...
            );
    end

    % Find common space/time nodes
    % Note: We assume the reference simulation has finer discretization both in
    % space and time; else we would have to check individually which simulation
    % has the finer discretization

    % Compute the number of skipped samples in time and space
    nSkipTime = simParsSim.h/simParsRef.h;
    nSkipSeg =  l_unit_Sim/l_unit_Ref;

    % get Reference simRes at common space time nodes
    simResRefRs.R = simDataRef.R(:,:,1:nSkipSeg:end,1:nSkipTime:end);
    simResRefRs.x = simDataRef.x(:,1:nSkipSeg:end,1:nSkipTime:end);
    simResRefRs.eta = simDataRef.eta(:,1:nSkipSeg:end,1:nSkipTime:end);
    simResRefRs.tout = simDataRef.tout(1:nSkipTime:end);


    %% Plot Preparation

    % Colors for plots with 3 and 9-dimensional data
    plotColors3 = lines(3);
    plotColors9 = lines(9);


    %% Get nodes / segments to plot
    % For beams with a high node number, don't plot the variables for
    % individual nodes/segments for all nodes, but only a maximum of n
    % segments / nodes. If there are more segments/nodes, skip some.


    nNodes = size(simDataSim.R, 3);
    nPlots = 3;

    if (nNodes-1) > nPlots
        plotNodes =  round(linspace(1, nNodes, nPlots));

        % Plot 2nd instead of first node since 1st node is always zero
        plotNodes(1) = 2;
    else
        plotNodes = 2:nNodes;
    end



    %% Plot configuration


    figHandles.comparisonConfig = figure( ...
        'Name', 'Comparison Configuration', 'NumberTitle','off');

    t = tiledlayout(2, length(plotNodes), 'TileSpacing','tight');
    title(t, 'Node Configuration $g_a$ (absolute)', 'interpreter', 'latex');

    ax = gobjects(2, length(plotNodes));

    for iTile = 1:length(plotNodes)

        iN = plotNodes(iTile);

        %%% Plot rotation matrices R
        ax(1, iTile) = nexttile(iTile);

        plot(simResRefRs.tout, reshape(simResRefRs.R(:,:,iN,:),[9, length(simResRefRs.tout)]), '-')
        hold on
        plot(simDataSim.tout, reshape(simDataSim.R(:,:,iN,:),[9, length(simDataSim.tout)]), '--')

        grid on
        xlim([simDataSim.tout(1), simDataSim.tout(end)]);
        title( sprintf('Rot. Matrix $R_{%d}$', iN ), 'Interpreter', 'latex')
        ax(1, iTile).ColorOrder = plotColors9;

        % legend(...
        %     [cellstr(num2str((1:9).', 'R_%d (Reference)'));...
        %     cellstr(num2str((1:9).', 'R_%d (Simulation'))],...
        %     'Interpreter', 'latex' ...
        %     );

        %%% Plot positions x
        ax(2, iTile) = nexttile(iTile + length(plotNodes));

        plot(simResRefRs.tout, reshape( simResRefRs.x(:, iN, :), 3, []), '-' );
        hold on
        plot(simDataSim.tout, reshape( simDataSim.x(:, iN, :), 3, []), '--' );

        grid on
        xlim([simDataSim.tout(1), simDataSim.tout(end)]);
        title( sprintf('Position $x_{%d}$', iN ), 'Interpreter', 'latex')
        legend('$x$', '$y$', '$z$', 'interpreter', 'latex');
        ax(2, iTile).ColorOrder = plotColors3;

        xlabel('time $t$ / s', 'interpreter', 'latex')

        legend(...
            [cellstr(num2str((1:3).', '$x_%d$ Ref.'));...
            cellstr(num2str((1:3).', '$x_%d$ Sim.'))],...
            'Interpreter', 'latex' ...
            );
    end

    % Synchronize axis limits of the subplot rows
    % linkaxes(ax(1,:));
    % linkaxes(ax(2,:));


    %% Plot velocities

    figHandles.comparisonEta = figure( ...
        'Name', 'Comparison Velocities', 'NumberTitle','off');

    t = tiledlayout(2, length(plotNodes), 'TileSpacing','tight');
    title(t, 'Node Velocities $\eta_a$ (body-fixed frame)', 'interpreter', 'latex');

    ax = gobjects(2, length(plotNodes));

    for iTile = 1:length(plotNodes)

        iN = plotNodes(iTile);

        %%% Plot rotational parts
        ax(1, iTile) = nexttile(iTile);

        plot(simResRefRs.tout, reshape( simResRefRs.eta(1:3, iN, :), 3, []), '-' );
        hold on;
        plot(simDataSim.tout, reshape( simDataSim.eta(1:3, iN, :), 3, []), '--' );

        grid on
        xlim([simDataSim.tout(1), simDataSim.tout(end)]);
        title( sprintf('Rot., Node %d', iN ), 'Interpreter', 'latex')
        legend(...
            [cellstr(num2str((1:3).', '$\\omega_%d$ Ref.'));...
            cellstr(num2str((1:3).', '$\\omega_%d$ Sim.'))],...
            'Interpreter', 'latex' ...
            );
        colororder(ax(1, iTile), plotColors3);


        %%% Plot translational parts
        ax(2, iTile) = nexttile(iTile + length(plotNodes));

        plot(simResRefRs.tout, reshape( simResRefRs.eta(4:6, iN, :), 3, []), '-' );
        hold on
        plot(simDataSim.tout, reshape( simDataSim.eta(4:6, iN, :), 3, []), '--' );

        grid on
        xlim([simDataSim.tout(1), simDataSim.tout(end)]);
        title( sprintf('Transl., Node %d', iN ), 'Interpreter', 'latex')

        legend(...
            [cellstr(num2str((1:3).', '$v_%d$ Ref.'));...
            cellstr(num2str((1:3).', '$v_%d$ Sim.'))],...
            'Interpreter', 'latex' ...
            );
        xlabel('time $t$ / s', 'interpreter', 'latex')
        colororder(ax(2, iTile), plotColors3);
    end

    % Synchronize axis limits of the subplot rows
    % linkaxes(ax(1,:));
    % linkaxes(ax(2,:));
end
