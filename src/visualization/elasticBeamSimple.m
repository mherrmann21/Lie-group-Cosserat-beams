classdef elasticBeamSimple < handle % Make it a handle function!
    % elasticBeam Class to draw an elastic beam in a 3D plot.
    %
    %
    % Maximilian Herrmann
    % Chair of Automatic Control
    % TUM School of Engineering and Design
    % Technical University of Munich
    %
    % Visualization of a geometrically exact beam
    % Direction of the beam length is in (local) Z-direction, i.e.,
    % the cross-sections lie in the XY planes of the local frames
    %
    % Properties and methods are mostly self-explanatory.


    properties

        % Nr. of segments
        nSeg = nan;

        % Array of SE3 node configuration matrices with dimension
        % (4, 4, nNodes)
        g = zeros(4,4,0);

        % Array of SE3 node configuration matrices with dimension
        % (4, 4, nNodesInterp),
        gInterp = zeros(4,4,0);

        % Beam cross-section width and height
        width
        height

        % Array with corners of cross sections in global coordinates
        % with dimensions (3, 4, nNodes)
        CSCorners = zeros(3, 4, 0);

        % Beam color
        col

        % Name of the beam (in the plot legend
        name

        % Interpolated segment length
        % Set to zero to disable spatial interpolation
        lInterp

        %%% Arrays with graphic objects

        % Center line plot3 object
        hCenterLine = gobjects(1);

        % Array with patches for cross-sections
        patchesCSec = gobjects(0);
    end

    methods
        function obj = elasticBeamSimple(g, options)
            %elasticBeam Construct an instance of this class
            %   ...

            arguments
                % Array of SE3 node configuration matrices with dimension
                % (4, 4, nNodes)
                g (4,4,:) double

                % Optional arguments

                % Beam width and height
                options.width  (1,1) {mustBeNumeric} = 0.015;
                options.height (1,1) {mustBeNumeric} = 0.015;

                options.color = lines(1);

                options.lInterp = 0.01;

                options.name = "";

                % Axis, in which to plot the beam
                % If none is given, use current axis
                options.axis = nan;
            end

            % Get nr of segments from data dimensions
            obj.nSeg = size(g,3) - 1;

            obj.g = g;
            obj.width   = options.width;
            obj.height  = options.height;
            obj.col     = options.color;
            obj.lInterp = options.lInterp;
            obj.name = options.name;

            if ~isgraphics(options.axis)
                options.axis = gca;
            end

            obj.drawBeam(options.axis);
        end

        function drawBeam(obj, ax)

            obj.computeCSCorners;
            if obj.lInterp
                obj.interpolateCenterLine;
            else
                obj.gInterp = obj.g;
            end

            for iN = 1:(obj.nSeg + 1)

                % Get absolute corner positions for current node
                cornersAbs = obj.CSCorners(:, :, iN);

                % Matrix with repeated first points (=first column) to close rectangle
                cornersAbsCl = [ cornersAbs, cornersAbs(:, 1)];

                % Plot cross-sections
                obj.patchesCSec(iN) = patch(ax, ...
                    cornersAbsCl(1,:), cornersAbsCl(2, :), cornersAbsCl(3, :), ...
                    obj.col, 'EdgeColor', obj.col, 'FaceAlpha', 0.3, ...
                    'HandleVisibility', 'off' ...
                    );
            end

            % Plot center line
            obj.hCenterLine = plot3(ax, ...
                squeeze(obj.gInterp(1,4,:)), ...
                squeeze(obj.gInterp(2,4,:)), ...
                squeeze(obj.gInterp(3,4,:)), ...
                'Color', obj.col ...
                );

            if ~isempty(obj.name)
                obj.hCenterLine.DisplayName = obj.name;
            end
        end

        function updateConfiguration(obj, g)
            % Update beam configuration
            arguments
                obj

                % Array of SE3 node configuration matrices with dimension
                % (4, 4, nNodes)
                g (4,4,:) double
            end

            obj.g = g;
            obj.computeCSCorners;
            if obj.lInterp
                obj.interpolateCenterLine;
            else
                obj.gInterp = obj.g;
            end

            for iNode = 1:(obj.nSeg + 1)

                % Get absolute corner positions for current node
                cornersAbs = obj.CSCorners(:, :, iNode);

                % Matrix with repeated first points (=first column) to close rectangle
                cornersAbsCl = [ cornersAbs, cornersAbs(:, 1)];

                % Update cross-sections
                obj.patchesCSec(iNode).XData = cornersAbsCl(1, :);
                obj.patchesCSec(iNode).YData = cornersAbsCl(2, :);
                obj.patchesCSec(iNode).ZData = cornersAbsCl(3, :);
            end

            % Update center line
            obj.hCenterLine.XData = obj.gInterp(1,4,:);
            obj.hCenterLine.YData = obj.gInterp(2,4,:);
            obj.hCenterLine.ZData = obj.gInterp(3,4,:);

        end

        function computeCSCorners(obj)
            % Computes the coordinates of the cross-section corners in
            % global frame

            % Matrix with corners in local frames
            % Vectors are in the rows; cross-sections are in the XY-plane
            % Order: Top left, top right, bottom right, bottom left
            corners = [
                -obj.width, +obj.height, 0
                -obj.width, -obj.height, 0
                +obj.width, -obj.height, 0
                +obj.width, +obj.height, 0
                ]*0.5;

            for iNode = 1:(obj.nSeg + 1)
                % Compute matrix with corner vectors in global frame
                % (Vectors in columns; makes computation easier)
                cornersAbs = obj.g(1:3,4,iNode) + obj.g(1:3,1:3,iNode) * corners';
                obj.CSCorners(:, :, iNode) = cornersAbs;
            end
        end

        function interpolateCenterLine(obj)
            % Spatial interpolation on SE3

            
            % Spatial query points (normalize beam length to 1)
            sInput = 0:1/obj.nSeg:1;            
            sQuery = 0:obj.lInterp:1;
            xi = computeDiscreteDeformations(obj.g)/(1/obj.nSeg);
            obj.gInterp = interpSimResSpaceSE3(obj.g,xi,sInput,sQuery);
        end
    end
end

