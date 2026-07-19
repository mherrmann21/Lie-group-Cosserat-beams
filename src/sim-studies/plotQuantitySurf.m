function figHandle = plotQuantitySurf(values, E, h, nSeg, modelNames, opts)
    arguments
        values      (:,:,:,:) double
        E           (:,1)
        h           (:,1) double
        nSeg        (:,1) double
        modelNames  (:,1) string

        opts
        % opts.figureName     (:,1) string
        % opts.figureTitle    (:,1) string
        % opts.axisLabel      (1,1) string
    end

    % Line styles and colors to be used by all plots
    plotMarkers = ["o", "*", "square", "x", "+"];
    lineStyles = ["-o", "-*", "-square", "-x", "-+", "-^", "-v"];
    lineColors = lines(length(lineStyles));

    % Make sure that the parameter ranges are large enough
    if length(nSeg) > 1 && length(h) > 1

        figHandle = figure('Name', opts.figureName);

        t = tiledlayout(1, length(E));
        title(t, opts.figureTitle, 'interpreter', 'latex');

        ax = gobjects(length(E));

        for iMat = 1:length(E)
            ax(iMat) = nexttile(t);

            for iSim = 1:length(modelNames)
                sh = surf(ax(iMat), ...
                    h, nSeg, squeeze(values(iMat, :,:, iSim))' );
                hold on

                sh.LineWidth = 1;
                sh.EdgeColor = lineColors(iSim, :);
                sh.FaceAlpha = 0.5;
                sh.FaceColor = lineColors(iSim, :);
                sh.FaceLighting = 'gouraud';
                sh.SpecularStrength = 0.5;
                sh.Marker    = plotMarkers(iSim);
                sh.MarkerSize = 10;
                sh.DisplayName = modelNames(iSim);
            end

            ax(iMat).XScale = 'log';
            ax(iMat).YScale = 'log';
            ax(iMat).ZScale = 'log';
            ax(iMat).YDir = 'reverse';

            light(ax(iMat), 'Position', [ ...
                0, ...
                0, ...
                max(ax(iMat).ZLim)*1.5 ...
                ]);

            title(sprintf('Mat=%d', iMat), 'Interpreter','latex');

            xlabel('time step $h$ / s', 'Interpreter','latex');
            ylabel('nr. of segments $n$', 'Interpreter','latex');
            zlabel(opts.axisLabel, 'Interpreter', 'latex');
            legend('Interpreter','latex', 'Location', 'best');
        end
        figHandle.NumberTitle = 'off';
    else
        % Return empty figure handle
        figHandle = gobjects(0);
    end
end