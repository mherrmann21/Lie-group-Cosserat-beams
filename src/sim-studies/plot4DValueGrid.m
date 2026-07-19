function [figHandle, ax] = plot4DValueGrid(values, dimInfo, figInfo, convLines)
    %% Plot 4D Data in a matrix of 2D subplots
    %
    % The data is plotted according to the dimensions:
    %  dimension 1: Y-Values in the individual subplots;
    %               dimInfo(1).xValues are the plot's x values;
    %               dimInfo(1).xAxisLabel is the plot x axis label
    %  dimension 2: Individual lines in each subplot;
    %               dimInfo(2).legendStrings(iLine) are the strings used in
    %               the plot legend;
    %               dimInfo(2).quantNameShort is used in () in the figure
    %               name
    %  dimension 3: Rows in the subplot matrix
    %               dimInfo(3).legendStrings(iRow) appears in the subplot
    %               title
    %  dimension 4: Columns in the subplot matrix
    %               dimInfo(4).legendStrings(iCol) appears in the subplot
    %               title

    % ToDo: Assert that the size of the data in the info struct
    % correspond to the size of values

    arguments
        values      (:,:,:,:)   double
        dimInfo     (:,1)       struct
        figInfo     (1,1)       struct

        % Struct with values for the convergence lines (in arrays);
        % dimensions are (Row, Col, nConvLines) for the subplot matrix
        convLines  (:,:,:)       struct
    end

    % Line styles and colors to be used by all plots
    lineStyles = repmat(["-o", "-*", "-square", "-x", "-+", "-^", "-v"], [1,ceil(size(values, 2)/7)]);
    lineColors = lines(length(lineStyles));

    % Initialize figure and subplot layout / axes
    figName = strcat( ...
        figInfo.quantNameShort, " / ", dimInfo(1).quantNameShort, ...
        " (", dimInfo(2).quantNameShort, ")" ...
        );

    figHandle = figure('Name', figName , 'NumberTitle', 'off');
    ax = gobjects(size(values, 3), size(values, 4));

    % Only plot if we don't have too many subplots
    if size(values, 4) < 15

        t = tiledlayout(size(values, 3), size(values, 4));

        title(t, ...
            strcat(figInfo.quantNameLong, " over ", dimInfo(1).quantNameLong), ...
            'interpreter', 'latex');

        plotEmpty = zeros(size(ax), 'logical');

        % Plot data
        for iRow = 1:size(values, 3)
            for iCol = 1:size(values, 4)

                ax(iRow, iCol) = nexttile(t);
                plotEmpty(iRow, iCol) = all(isnan(values(:, :, iRow, iCol)), 'all');

                if ~plotEmpty(iRow, iCol)
                    for iLine = 1:size(values, 2)

                        if ~all(isnan(values(:, iLine, iRow, iCol)))

                            plot(ax(iRow, iCol), ...
                                dimInfo(1).xValues, values(:, iLine, iRow, iCol), ...
                                lineStyles(iLine), ...
                                'Color', lineColors(iLine, :), 'MarkerSize', 9, ...
                                'LineWidth', 1, ...
                                ...%'MarkerFaceColor', opts.plotColors(iSim, :), 'MarkerEdgeColor', opts.plotMarkerEdgeColors{iSim},...
                                'DisplayName', dimInfo(2).legendStrings(iLine) ...
                                );
                            hold on;

                        end
                    end

                    % Plot convergence bound lines
                    % ToDo: Could be extended to plot multiple conv lines in
                    % one plot
                    if ~isempty(convLines)
                        for iConv = 1:size(convLines.a,3)

                            % Get coefficients for all integrators
                            a = convLines.a(iRow, iCol, iConv);
                            b = convLines.b(iRow, iCol, iConv);

                            if ~isnan(a) && ~isnan(b)
                                convVals = a*(dimInfo(1).xValues).^b;

                                equString = sprintf('$%.3g \\, x^{%.3g}$', a, b);

                                plot(ax(iRow, iCol), ...
                                    dimInfo(1).xValues, convVals, ...
                                    '--', ...
                                    'Color', lineColors(iConv, :), ...
                                    'DisplayName', equString ...
                                    );
                            end
                        end
                    end
                    legend('Interpreter','latex', 'Location', 'best');
                end
                ax(iRow, iCol).XScale = 'log';
                ax(iRow, iCol).YScale = 'log';
                if numel(dimInfo(1).xValues) > 1
                    ax(iRow, iCol).XLim = sort(dimInfo(1).xValues([1,end]));
                end
                %ax(iRow, iCol).DataAspectRatio = [1,1,1];

                ax(iRow, iCol).TickLabelInterpreter = "latex";

                xlabel(dimInfo(1).xAxisLabel, 'Interpreter','latex');
                ylabel(figInfo.yAxisLabel, 'Interpreter', 'latex');
                title( ...
                    strcat(dimInfo(3).legendStrings(iRow), ", ", dimInfo(4).legendStrings(iCol)), ...
                    'Interpreter','latex' ...
                    );
                grid on;

            end
        end
        % Synchronize axis limits
        % Note: Specifically only synchronize non-empty plots, since these
        % would corrupt the correct axis limits (workaround)
        if ~isempty(ax(~plotEmpty))
            linkaxes(ax(~plotEmpty));
        end
    end
end