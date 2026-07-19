function [saveDirCase,saveDirAll] = prepareSimStudyFolderStructure(opts,beamPars,h,nSeg)
    %% Prepare the output folder structure for a simulation study
    % 
    % Inputs: See arguments block
    % Outputs: String array with file paths to the individual simulation
    % case folders

    arguments
        % Settings struct with fields
        %    folderSuffix
        %    resultsDir
        opts        (1,1) struct

        % Vectors with the parameter values for the simulation study
        beamPars    (:,1) 
        h           (:,1) double
        nSeg        (:,1) double
    end


    disp('Preparing Output Folder Structure...')

    % Get subfolder for all simulation results
    subFolderAll = strcat( string(datetime, 'yyMMdd_HHmm'), opts.folderSuffix);

    % Full path for all simulation results
    saveDirAll = fullfile(opts.resultsDir, subFolderAll);

    saveDirCase = strings( ...
        length(beamPars), ...
        length(h), ...
        length(nSeg));

    for iMat = 1:length(beamPars)
        for ih = 1:length(h)
            for iSeg = 1:length(nSeg)
                % Get subfolder for current simulation case
                subFolderCase = sprintf( ...
                    'Case %d-%d-%d - Mat=%d h=%.2E nSeg=%d', ...
                    iMat, ih, iSeg, iMat, h(ih), nSeg(iSeg));

                % Get path for the results directory
                saveDirCase(iMat, ih, iSeg) = ...
                    fullfile(saveDirAll, subFolderCase);

                % Create output folder
                if ~isfolder( saveDirCase(iMat, ih, iSeg) )
                    mkdir( saveDirCase(iMat, ih, iSeg) );
                end
            end
        end
    end
end