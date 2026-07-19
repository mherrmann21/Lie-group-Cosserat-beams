classdef beamParamsDiscrete
    % Parameter class to represent the properties of a discretized beam,
    % i.e., discrete inertial / stiffness parameters
    %
    % Maximilian Herrmann
    % Chair of Automatic Control
    % TUM School of Engineering and Design
    % Technical University of Munich


    properties

        nAllwd    (1,1)   double   % Number of allowed deformation modes
        nSeg      (1,1)   double   % Nr. of Segments
        nNodes    (1,1)   double   % Nr. of Nodes (nNodes = nSeg + 1)
        l         (1,1)   double   % Segment Length

        % (Reduced) stiffness matrix with dimensions (nAllwd,nAllwd)
        CgenRed   (:,:)   double

        MgenNode  (6,6,:) double   % Generalized inertia tensors for the discrete nodes
        mNode     (1,:)   double   % Total mass of the discrete nodes

        % Discrete dissipation matrix of a single segment for linear
        % dissipation; dimensions are (nAllwd, 6), i.e., it is computed
        % with Ba.' * diag(dLin) * Ba
        DLinRed   (:,:)   double

        % Pre-computed dissipation matrix (linear dissipation) for the
        % entire system; dimensions are (nAllwd*nSeg, nAllwd*nSeg)
        DLinRedSys (:,:) double

        % Vector of discrete dissipation coefficients for quadratic dissipation
        dQuad     (6,1)   double

        % Vector of body masses (of attached bodies); dimensions (nNodes,1)
        m_a         (:,1)   double

        % Array of position vectors from the cross-section frame to the COM
        % frame of attached rigid bodies
        x_a         (3,:) double
    end

    methods

        function obj = beamParamsDiscrete(simPars,beamPars,Ba)
            if nargin == 0
                % Return empty object
            elseif nargin == 3
                % Compute values from given arguments
                obj = obj.computeDiscreteParams(simPars,beamPars,Ba);
            else
                error('Argument Count not Supported.')
            end

        end

        function obj = computeDiscreteParams(obj,simPars,beamPars,Ba)

            % Nr. of nodes and segments
            obj.nNodes = size(simPars.g0, 3);
            obj.nSeg   = obj.nNodes - 1;
            obj.nAllwd = size(Ba, 2);

            % Segment length
            obj.l = beamPars.L / obj.nSeg;

            % Reduced stiffness matrix
            obj.CgenRed  = obj.l * Ba.' * beamPars.Cgen * Ba;

            % Dissipation variables
            obj.DLinRed    = obj.l * Ba.' * diag(beamPars.d) * Ba;
            obj.DLinRedSys = obj.l * diag(repmat(Ba.' * beamPars.d, [obj.nSeg,1]));
            obj.dQuad      = obj.l * beamPars.dq;

            %%% Compute generalized inertia tensors and mass values for the
            % discrete nodes

            % Check that both inertia tensor and mass arrays are either
            % both empty or have the correct dimensions
            assert( size(simPars.M_a, 3) == size(simPars.m_a, 1), ...
                ['Inconsistent array sizes for masses and inertia tensors of attached bodies. ' ...
                'To include attached bodies, both arrays must have values and the correct dimensions.'] ...
                );          
            assert( size(simPars.M_a, 3) == size(simPars.x_a, 2), ...
                ['Inconsistent array sizes for position vectors and inertia tensors of attached bodies. ' ...
                'To include attached bodies, both arrays must have values and the correct dimensions.'] ...
                );

            if isempty(simPars.M_a)
                obj.MgenNode = zeros(6,6,obj.nNodes);
                obj.mNode    = zeros(1,obj.nNodes);
                obj.m_a      = zeros(1,obj.nNodes);
                obj.x_a      = zeros(3,obj.nNodes);
            else
                obj.MgenNode = simPars.M_a;
                obj.mNode    = simPars.m_a;
                obj.m_a      = simPars.m_a;                
                obj.x_a      = simPars.x_a;                
            end

            obj.MgenNode = obj.MgenNode + obj.l * cat(3, ...
                0.5*beamPars.Mgen, ...
                repmat(beamPars.Mgen, 1, 1, obj.nNodes-2), ...
                0.5*beamPars.Mgen ...
                );

            obj.mNode = obj.mNode ...
                + obj.l * [0.5, ones(1,obj.nNodes-2), 0.5] * beamPars.m;

        end

    end

end
