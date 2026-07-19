classdef beamSimulation
    % Class to store everything needed to fully simulate a given beam model
    %
    % Maximilian Herrmann
    % Chair of Automatic Control
    % TUM School of Engineering and Design
    % Technical University of Munich

    properties
        simModel    (1,1) beamSimModel

        simRes      (1,1) beamSimRes

        simPars     (1,1) beamSimPars

        beamPars    (1,1) beamParams
    end

    methods
        function obj = simulateModel(obj, options)
            % Simulate the model using the stored parameters
            % The function handle must call the corresponding simulation
            % function that accepts the inputs as defined in the function
            % call below
            arguments
                obj

                % If true, integration time is measured with timeit. This
                % requires additional simulation runs but is more accurate.
                options.accurateTiming (1,1) logical = false;

                % Whether or not to print simulation metadata to the
                % console
                options.consoleOutput (1,1) logical = true;
            end

            % ToDo: Validate that all necessary parameters are available

            %%% Check if the dimensions of the output arrays for desired simulation
            % and discretization parameters are small enough for C++ Code,
            % see https://www.mathworks.com/help/coder/ug/array-size-restrictions-for-code-generation.html
            % We check for the 3x3 rotation matrix array since this is the
            % largest one.
            % Nr. of integration steps
            nSteps = round( obj.simPars.tEnd / obj.simPars.h );
            nNodes = size(obj.simPars.g0, 3);
            if (3*3*nSteps*nNodes) > double(intmax())
                % Dimensions too large: Skip simulation
                obj.simRes.metaDataSim.exitCode = 0;
                fprintf('Integration skipped (%s): Output array dimensions too large.\n', obj.simModel.modelName);
            else
                % Dimensions okay: Run integration
                if options.consoleOutput
                    fprintf('Starting integration (%s)...\n', obj.simModel.modelName);
                end

                % Execute simulation, measure time
                tic
                [obj.simRes.simData, obj.simRes.metaDataSteps] = obj.simModel.funHandle(obj.beamPars, obj.simPars, obj.simModel);
                obj.simRes.metaDataSim.TotalTime = toc;

                % Get exit code / check if simulation was successful
                obj.simRes.metaDataSim.exitCode = getSimulationExitCode( ...
                    obj.simRes.simData, obj.simRes.metaDataSteps, ...
                    obj.simPars, obj.simModel.solverConfig ...
                    );

                % Run simulation again with timing (if desired) for
                % accurateTiming times
                if obj.simRes.metaDataSim.exitCode && options.accurateTiming
                    timingFun = @() obj.simModel.funHandle(obj.beamPars, obj.simPars, obj.simModel);
                    obj.simRes.metaDataSim.TotalTime = timeit(timingFun);
                end

                obj.simRes.metaDataSim.StepTimeMean = ...
                    obj.simRes.metaDataSim.TotalTime / length(obj.simRes.simData.tout);

                % Compute some metadata values of the simulation
                obj.simRes = obj.simRes.getSimMetaData();

                % Print some information
                if options.consoleOutput
                    if obj.simRes.metaDataSim.exitCode && options.accurateTiming
                        fprintf('   Total integration time (timeit):       %f s\n', obj.simRes.metaDataSim.TotalTime);
                    else
                        fprintf('   Total integration time (tictoc):       %f s\n', obj.simRes.metaDataSim.TotalTime);
                    end
                    fprintf('   Nr. of successful time steps:          %d\n', length(obj.simRes.simData.tout));
                    fprintf('   Simulation end time:                   %f s\n', obj.simRes.simData.tout(end));
                    fprintf('   Time step:                             %f ms = 2^%.2f s\n', obj.simPars.h*1000, log2(obj.simPars.h));
                    fprintf('   Approx. comp. time per step:           %f ms\n', ...
                        obj.simRes.metaDataSim.TotalTime / length(obj.simRes.simData.tout)*1e3);
                    fprintf('   Iteration Count (Avg/Min/Max):         %.4f / %2d / %2d\n', ...
                        obj.simRes.metaDataSim.ImplicitIterations.mean, ...
                        obj.simRes.metaDataSim.ImplicitIterations.min, ...
                        obj.simRes.metaDataSim.ImplicitIterations.max );
                    fprintf('   Residual (Avg/Min/Max):                %.4g / %.4g / %.4g\n', ...
                        obj.simRes.metaDataSim.ImplicitError.mean, ...
                        obj.simRes.metaDataSim.ImplicitError.min, ...
                        obj.simRes.metaDataSim.ImplicitError.max );
                    fprintf('   Nr. of target solver error violations: %d (%.2f%% of overall steps)\n', ...
                        sum(obj.simRes.metaDataSteps.ExitFlag == 1), sum(obj.simRes.metaDataSteps.ExitFlag == 1) / length(obj.simRes.simData.tout) *100);
                    fprintf('   Simulation exit code:                  %d\n', obj.simRes.metaDataSim.exitCode);
                end
            end
        end

        function obj = computeEnergyEvolution(obj)
            % Compute the energy evolution of the simulation results

            arguments
                obj
            end
            % First compute discrete deformations in case they aren't
            % present yet
            obj.simRes.computeDeformations(obj.beamPars);

            % Then compute the actual energy evolution
            obj.simRes.computeEnergyEvolution(obj.simPars.xiRef, obj.beamPars, obj.simPars);
        end

        function obj = computeAllResultsData(obj)
            % Compute all missing simulation results quantities
            obj.simRes.computeAllResultsData(obj.simPars.xiRef, obj.beamPars, obj.simPars);
        end

        function obj = clearSimData(obj)
            % Deletes all simulation results data except for the overall
            % simulation metadata;
            % this is done to be able to save the object without the large
            % simulation data

            obj.simRes.simData = beamSimData;
            obj.simRes.metaDataSteps  = beamSimMetaDataSteps;
            obj.simRes.E = beamSimEnergies;
        end


        function saveObject(obj, filepath)
            % Save this object (the variable) to a file
            arguments
                obj

                % Full path to the target file (including file name and
                % file ending)
                filepath (1,1) string
            end

            % Use .mat file version 7.3 to be able to save larger
            % variables
            try
                save(filepath, 'obj', '-mat', '-v7.3');
            catch ME
                warning( ...
                    ME.identifier, ...
                    'Could not save object:\n %s', ME.message ...
                    );
            end
        end


        function figHandles = plotSimResults(obj, options)
            % Plot the simulation results of this model / simulation and
            % (optionally) save them to file
            arguments
                obj                   (1,1) beamSimulation

                % Time step at which the results are plotted
                % IMPORTANT: This resamples the results in simRes at the
                % given sample time; i.e. the fine result data is LOST!
                options.hPlot         (1,1) double = nan;

                % Compute energy for plotting
                options.computeEnergy (1,1) logical = false;

                options.saveDir       (1,1) string = ""
                options.savePlotsJPEG (1,1) logical = false
                options.savePlotsFig  (1,1) logical = false
            end

            % Interpolate results, if required
            if ~isnan(options.hPlot)
                obj.simRes.interpolateSimResTime(options.hPlot);
            end

            % Compute energy, if required
            if options.computeEnergy
                obj.simRes.computeEnergyEvolution(obj.simPars.xiRef, obj.beamPars, obj.simPars);
            end

            % Plot results
            figHandles = plotBeamSimRes( ...
                obj.simRes, obj.simPars.xiRef, obj.simModel.modelName ...
                );

            % Plot reference and initial condition
            figHandles.RefConfig     = visualizeBeamConfig(obj.simPars.gRef, obj.beamPars, 'Reference Configuration');
            figHandles.InitialConfig = visualizeBeamConfig(obj.simPars.g0,   obj.beamPars, 'Initial Configuration');

            % Plot external forces
            %figHandles.ExtForces = plotExtForces(obj.simPars);

            drawnow;

            % Save to file if required
            if options.savePlotsFig || options.savePlotsJPEG
                disp('Saving Simulation Results Plots...')

                % Get model name without spaces
                namePrefix = replace(obj.simModel.modelName, " ", "_");

                saveFiguresFromStruct( ...
                    figHandles, options.saveDir, ...
                    "namePrefix", strcat(namePrefix, " "), ...
                    "saveFig",  options.savePlotsFig, ...
                    "saveJPEG", options.savePlotsJPEG ...
                    );
            end
        end


        function figHandles = plotReferenceComparison(obj, beamMdlRef, options)
            % Plot the comparison between the simulation results and a
            % reference model and (optionally) save them to file
            arguments
                obj                   (1,1) beamSimulation

                % beamSimModel object with the reference simulation results
                beamMdlRef            (1,1) beamSimulation

                options.saveDir       (1,1) string = ""
                options.savePlotsJPEG (1,1) logical = false
                options.savePlotsFig  (1,1) logical = false
            end

            % Plot comparison
            figHandles = plotBeamReferenceComparison( ...
                obj.simPars, obj.simRes.simData, ...
                beamMdlRef.simPars, beamMdlRef.simRes.simData ...
                );

            % Save to file if required
            if options.savePlotsFig || options.savePlotsJPEG
                disp('Saving Reference Comparison Plots...')

                % Get model name without spaces
                namePrefix = replace(obj.simModel.modelName, " ", "_");

                saveFiguresFromStruct( ...
                    figHandles, options.saveDir, ...
                    "namePrefix", strcat(namePrefix, " "), ...
                    "saveFig",  options.savePlotsFig, ...
                    "saveJPEG", options.savePlotsJPEG ...
                    );
            end
        end
    end
end
