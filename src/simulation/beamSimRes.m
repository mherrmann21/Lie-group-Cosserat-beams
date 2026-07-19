classdef beamSimRes
    % Class to store all simulation results of a beam simulation
    %
    % Maximilian Herrmann
    % Chair of Automatic Control
    % TUM School of Engineering and Design
    % Technical University of Munich

    properties
        simData       (1,1) beamSimData
        metaDataSteps (1,1) beamSimMetaDataSteps
        metaDataSim   (1,1) beamSimMetaDataSim
        E             (1,1) beamSimEnergies
    end

    methods

        function obj = computeAllResultsData(obj, xiRef, beamPars, simPars)
            % Compute all missing quantities

            arguments
                obj

                % Array of reference deformations with dimensions (6, nSeg)
                xiRef   (6,:)   double

                % Standard beam parameter struct
                beamPars  (1,1)   beamParams

                simPars   (1,1)   beamSimPars
            end

            if isempty(obj.simData.xi)
                obj = obj.computeDeformations(beamPars);
            end
            if isempty(obj.simData.xDot)
                obj = obj.computeVelocities();
            end
            if isempty(obj.metaDataSteps.orthErrorR)
                obj = obj.computeOrthError();
            end
            if isempty(obj.E.H)
                obj = obj.computeEnergyEvolution(xiRef, beamPars, simPars);
            end
        end


        function obj = computeDeformations(obj, params)
            % Compute discrete deformations / strains + stresses
            % Note: The computed values for xi correspond to discrete updates, not
            % gradients; thus, they are divided by segment length l afterwards
            % to get the gradients

            arguments
                obj
                params (1,1) beamParams
            end

            nSeg = size(obj.simData.R, 3) - 1;

            disp('   Computing Deformations...')
            obj.simData.xi = simResComputeDeformations_mex( ...
                obj.simData.tout, obj.simData.R, obj.simData.x, obj.simData.eta ...
                ) / (params.L / nSeg) ;
        end


        function obj = computeVelocities(obj)
            % Compute velocities in inertial frame
            disp('   Computing Velocities...')
            obj.simData.xDot = simResComputeVelocities_mex( ...
                obj.simData.tout, obj.simData.R, obj.simData.x, obj.simData.eta ...
                );
        end


        function obj = computeOrthError(obj)
            % Compute orthogonality error of rotation matrices

            disp('   Computing Orthogonality Error...')
            obj.metaDataSteps.orthErrorR = simResComputeOrthError_mex(obj.simData.R);
        end


        function obj = computeEnergyEvolution(obj, xiRef, beamPars, simPars)
            % Compute the energy evolution of the simulation results

            arguments
                obj

                % Array of reference deformations with dimensions (6, nSeg)
                xiRef   (6,:)   double

                % Standard beam parameter struct
                beamPars  (1,1)   beamParams

                simPars   (1,1)   beamSimPars
            end

            disp('   Computing Energy Evolution...')
            obj.E = computeBeamEnergyEvolution_mex( ...
                obj.simData.x, obj.simData.eta, obj.simData.xi, ...
                xiRef, beamPars, simPars ...
                );
        end

        function obj = getSimMetaData(obj)
            % Compute metadata for the entire simulation from the metadata
            % of the individual steps and simulation results
            % data

            obj.metaDataSim.ImplicitIterations.min  = min( obj.metaDataSteps.ImplicitIterations(:));
            obj.metaDataSim.ImplicitIterations.max  = max( obj.metaDataSteps.ImplicitIterations(:));
            obj.metaDataSim.ImplicitIterations.mean = mean(obj.metaDataSteps.ImplicitIterations(:), 'omitnan');

            obj.metaDataSim.ImplicitError.min       = min(  abs( obj.metaDataSteps.ImplicitError(:) ));
            obj.metaDataSim.ImplicitError.max       = max(  abs( obj.metaDataSteps.ImplicitError(:) ));
            obj.metaDataSim.ImplicitError.mean      = mean( abs( obj.metaDataSteps.ImplicitError(:) ), 'omitnan');

            obj.metaDataSim.TotalIterations         = sum( obj.metaDataSteps.ImplicitIterations(:), 'omitmissing');
        end

        function obj = interpolateSimResTime(obj, hResample)

            if hResample < obj.simData.tout(end)

                tQuery = (obj.simData.tout(1):hResample:obj.simData.tout(end))';

                nNodes = size(obj.simData.R,3);
                nSeg   = nNodes-1;

                % Use linear interpolation in time
                obj.simData.R = interpn( ...
                    1:3, 1:3, 1:nNodes, obj.simData.tout, ...
                    obj.simData.R, ...
                    1:3, 1:3, 1:nNodes, tQuery);

                obj.simData.x = interpn( ...
                    1:3, 1:nNodes, obj.simData.tout, ...
                    obj.simData.x, ...
                    1:3, 1:nNodes, tQuery ...
                    );

                obj.simData.eta = interpn( ...
                    1:6, 1:nNodes, obj.simData.tout, ...
                    obj.simData.eta, ...
                    1:6, 1:nNodes, tQuery ...
                    );

                obj.simData.xi = interpn( ...
                    1:6, 1:nSeg, obj.simData.tout, ...
                    obj.simData.xi, ...
                    1:6, 1:nSeg, tQuery ...
                    );

                if ~isempty(obj.metaDataSteps.orthErrorR)
                    obj.metaDataSteps.orthErrorR = interpn( ...
                        1:size(obj.metaDataSteps.orthErrorR, 1), obj.simData.tout, ...
                        obj.metaDataSteps.orthErrorR, ...
                        1:size(obj.metaDataSteps.orthErrorR, 1), tQuery ...
                        );
                end

                % 2D-Arrays don't need interpn; we just have to take care of
                % the transposes
                if ~isempty(obj.metaDataSteps.ImplicitIterations)
                    obj.metaDataSteps.ImplicitIterations = interp1( ...
                        obj.simData.tout, ...
                        obj.metaDataSteps.ImplicitIterations.', ...
                        tQuery ...
                        ).';
                end

                if numel(obj.metaDataSteps.ImplicitError) > 1
                    obj.metaDataSteps.ImplicitError = interp1( ...
                        obj.simData.tout, ...
                        obj.metaDataSteps.ImplicitError.', ...
                        tQuery ...
                        ).';
                end

                if ~isempty(obj.metaDataSteps.ExitFlag)
                    obj.metaDataSteps.ExitFlag = interpn( ...
                        obj.simData.tout, ...
                        obj.metaDataSteps.ExitFlag.', ...
                        tQuery ...
                        ).';
                end
                obj.simData.tout = tQuery;
            end
        end
    end

end
