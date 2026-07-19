function gQuery = interpSimResSpaceSE3(gInput,xiInput,sInput,sQuery)
    %% Geometrically correct spatial interpolation of beam configurations on SE3
    % Inputs: See below.
    %
    % Maximilian Herrmann
    % Chair of Automatic Control
    % TUM School of Engineering and Design
    % Technical University of Munich

    arguments
        % Array of SE3 configuration matrices of the simulation results
        % with dimensions (4,4,nNodes,nSteps)
        gInput  (4,4,:,:) double

        % Array of discrete deformations with dimensions (6,nSeg,nSteps)
        xiInput (6,:,:) double

        % Vector of the node positions s along the beam (for input
        % discretization)
        sInput  (:,1) double

        % Vector of desired query points (node positions)
        sQuery  (:,1) double
    end

    % Output array
    gQuery = zeros(4,4,length(sQuery), size(gInput,4));

    for iStep = 1:size(gInput,4)
        for iSQuery = 1:length(sQuery)

            % Check whether we are the end length
            if sQuery(iSQuery) == sInput(end)
                % End reached:
                % Simply take last node configuration
                gQuery(:,:,iSQuery,iStep) = gInput(:,:,end,iStep);
            else
                % We are somewhere on the beam (including s=0):
                % Compute via interpolation

                % Get segment in the original discretization, to which the current s
                % belongs:
                % Get index of first node (node position) that is larger than the
                % current position
                iNextNode = find(sInput>sQuery(iSQuery),1);

                % Configuration of the previous node (at the beginning of the segment)
                g0 = gInput(:,:,iNextNode-1,iStep);

                % Get position of the current query point w.r.t. the previous node of
                % the initial discretization
                lRel = sQuery(iSQuery) - sInput(iNextNode-1);

                % Compute configuration
                gQuery(:,:,iSQuery,iStep) = g0*caySE3(lRel * xiInput(:,iNextNode-1, iStep));
            end
        end
    end
end
