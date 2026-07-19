function saveFiguresFromStruct(figHandles, saveDir, options)
    %% saveFiguresFromStruct
    % Saves figures to file.
    % The figure handles are given as fields of a struct (array).
    % The name of the file is the corresponding struct field name;
    % optionally, a prefix can be added before the struct field name.
    %
    % ToDo: Add options to the saved figures, if needed
    %
    % Note: For plots with large amounts of data (lots of data points),
    % the generated files might be excessively large with default settings.
    % Hence, make sure that the default .mat file format in the MATLAB 
    % preferences is 7.3. (Preferences / General / Mat Files)
    % https://www.mathworks.com/matlabcentral/answers/1575343-error-using-savefig-and-saveas
    %
    % Maximilian Herrmann
    % Chair of Automatic Control
    % TUM School of Engineering and Design
    % Technical University of Munich

    arguments
        % Struct containing the figure handles as fields
        figHandles (1,1) struct

        % String with full path to the folder, where the files should be saved
        saveDir    (1,1) string

        % Optional prefix for the file name, which is added before the
        % struct field name
        options.namePrefix (1,1) string = "";

        options.saveJPEG (1,1) logical = true % Save figure as jpeg?
        options.saveFig  (1,1) logical = true % Save figure as fig?
    end

    plotFields = fieldnames(figHandles);

    for iPlot = 1:length(plotFields)
        if ~isempty( figHandles.(plotFields{iPlot}) )

            % Get full file path
            filePath = fullfile( ...
                saveDir, strcat(options.namePrefix, plotFields{iPlot}));

            try
                % Save as figure
                if options.saveFig
                    savefig( figHandles.(plotFields{iPlot}), ...
                        sprintf('%s.fig', filePath) ...
                        );
                end

                % Save as jpeg
                if options.saveJPEG
                    exportgraphics( ...
                        figHandles.(plotFields{iPlot}), ...
                        sprintf('%s.jpeg', filePath));
                end

            catch ME
                warning( ...
                    ME.identifier, ...
                    'Could not save figure:\n %s', ...
                    ME.message ...
                    );
            end
        end
    end
end
