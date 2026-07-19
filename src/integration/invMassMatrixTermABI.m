function q_ddot = invMassMatrixTermABI( ...
        SNode, g_ij, x, CayRTDMassMatrix_k, h ...
        )
    %% Compute the inverse mass matrix multiplied with a given vector
    % I.e., the result of the product (M^-1)*x, where M is the system's
    % mass matrix and x is a given vector.
    % The system is a Serial-Kinematic Multibody System.
    % This term is computed with the Articulated Body Algorithm from
    % [Kim12], table 3, and substituting zero velocities as described in
    % [Lee+20], sec. 3.2.
    %
    %  * Implementation is done with the following indexing convention
    %    (same as in [Kim12]):
    %    Joint i = joint connecting body i with its parent i-1
    %  * Here, we consider only systems with serial structure, i.e., every
    %    body has only one child body attached to it
    %  * The following quantities in the original algorithm are set to
    %    zero / omitted:
    %      * q_dot (generalized velocities)
    %      * eta   (absolute velocities in se3)
    %      * F_ext (external forces in se*3)
    %
    %  * The input / output is not in vectors, but in arrays of the form
    %    (nAllwd, nSeg)
    %
    % Maximilian Herrmann
    % Chair of Automatic Control
    % TUM School of Engineering and Design
    % Technical University of Munich

    arguments (Input)
        % Array of joint Jacobians with size (6, nAllwd, nSeg)
        SNode   (6,:,:) double

        % Array of relative node configurations
        % Index a corresponds to the update from the parent node a-1 to
        % the current node a
        g_ij    (4,4,:) double

        % Vector x, which is multiplied with the inverse mass matrix;
        % given as an array with size (nAllwd, nSeg)
        % In the original algorithm, these are the joint torques tau
        x       (:, :)  double

        % Array of discrete velocities from last time step with size
        % (6,nSeg) (NOTE: Excludes the first, fixed node!)
        % (needed for the node inertia matrices for the integrators)
        CayRTDMassMatrix_k   (6,6,:) double

        % Time step
        h        (1,1) double

    end

    arguments (Output)
        % Array of (relative) joint accelerations; same size as x
        % These are the desired values of the vector (M^-1)*x
        q_ddot  (:,:)   double
    end


    %% Get system properties

    % Nr. of DoFs per joint
    nAllwd = size(SNode,2);

    % Nr. of bodies / nodes
    nSeg   = size(SNode,3);


    %% Pass 1 (Forward 0 -> n): Configuration, aux variables, ...
    % Skipped for q_dot = 0 -- Nothing to do here


    %% Pass 2 (Backward n -> 0): Compute I_hat and B_hat
    % Note: This pass is implemented for the general case; i.e. arbitrary
    % q_dot!

    % Array of Articulated Body Inertias
    I_hat = zeros(6,6,nSeg);

    % Array of Bias Forces
    B_hat = zeros(6,nSeg);

    % Array for Intermediate variable Psi
    Psi = zeros(nAllwd, nAllwd, nSeg);

    % Array to store values for AdInv operator (to save computation time)
    AdInv_g_ij = zeros(6,6,nSeg);

    % Actual loop
    for ii = nSeg:(-1):1

        % Fixed term of I_hat -- only the generalized inertia matrix
        % (Specifically computed here for the variational integrators)
        I_hat(:,:,ii) = CayRTDMassMatrix_k(:,:,ii)/h;

        % B_hat has no fixed term considering zero velocities and F_ext = 0

        % Compute Ad operator for node
        AdInv_g_ij(:,:,ii) = lAdSE3Inv(g_ij(:,:, ii));

        % Add terms from child body
        % Note: We only consider one child body here (which means we don't
        % need a loop); if we want to consider multiple child bodies, the
        % computation must be done in a loop for each body
        if ii < nSeg
            ik = ii + 1;

            % Intermediate variable Pi ([Kim12], Table 3)
            Pi_k = ...
                + I_hat(:,:,ik) ...
                - I_hat(:,:,ik) * SNode(:,:,ik) * Psi(:,:,ik) * SNode(:,:,ik).' * I_hat(:,:,ik);

            % Intermediate variable beta ([Kim12], Table 3)
            beta_k = ...
                + B_hat(:,ik) ...
                + I_hat(:,:,ik) * SNode(:,:,ik) * Psi(:,:,ik) ...
                * ( x(:,ik) - SNode(:,:,ik).' *B_hat(:,ik) );

            % Add terms
            I_hat(:,:,ii) = ...
                + I_hat(:,:,ii) ...
                + AdInv_g_ij(:,:,ik).' * Pi_k * AdInv_g_ij(:,:,ik);

            B_hat(:,ii)   = ...
                B_hat(:,ii) ...
                + AdInv_g_ij(:,:,ik).' * beta_k;
        end

        % Intermediate variable Psi
        Psi(:,:,ii) = inv( SNode(:,:,ii).' * I_hat(:,:,ii) * SNode(:,:,ii) );
    end


    %% Pass 3 (Forward 0 --> n): Compute Accelerations

    % Relative joint accelerations
    q_ddot  = zeros(nAllwd, nSeg);

    % Absolute acceleration of the parent body in the chain
    eta_dot_i0 = zeros(6,1);

    for ii = 1:nSeg

        % Relative acceleration
        q_ddot(:,ii) = ...
            + Psi(:,:,ii) * (...
            + x(:,ii) ...
            - SNode(:,:,ii).' * I_hat(:,:,ii) ...
            * AdInv_g_ij(:,:,ii) * eta_dot_i0 ...
            - SNode(:,:,ii).' * B_hat(:,ii) ...
            );

        % Absolute acceleration of current body i
        % (Directly assigned to variabl for parent body for the next
        % iteration)
        eta_dot_i0 = ...
            + AdInv_g_ij(:,:,ii) * eta_dot_i0 ...
            + SNode(:,:,ii) * q_ddot(:,ii);
    end

end
