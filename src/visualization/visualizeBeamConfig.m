function [figHandle, beamVisHandle] = visualizeBeamConfig(g, params, Name, options)
    %% Visualize a beam configuration and its deformation gradients
    % Returns the figure and elastic-beam graphics handles.
    arguments
        % Array of SE3 node configuration matrices with dimension 
        % (4, 4, nNodes)
        g (4,4, :) double

        % Struct with the beam parameters; needed for beam length and width
        params (1,1) beamParams

        % Name of the configuration (character array or string)
        Name (:,:)
        
        options.showFrames (1,1) logical = false
        options.showLabels (1,1) logical = false
    end

    nNodes = size(g, 3);
    nSeg = nNodes - 1;
    l = params.L / nSeg;

    beamW = params.geom.W;
    beamH = params.geom.H;

    % Compute discrete deformations
    % Note: Round slightly to remove numeric errors (just for
    % visualization!)
    errorLevel = 1e10;
    xi = round(computeDiscreteDeformations(g)/l*errorLevel)/errorLevel;

    figHandle = figure('Name', Name);
    tiledlayout(2,2);
    nexttile(1,[2,1]);

    init3Dplot('createFigure', false);
    xlabel('x'); ylabel('y'); zlabel('z');
    axis equal

    beamVisHandle = elasticBeam( ...
        g, 'showLabels', options.showLabels, 'showFrames', options.showFrames, ...
        'width', beamW, 'height', beamH ...
        );
    title(Name, 'Interpreter','latex');

    % Plot discrete deformations (rotational)
    ax = nexttile();

    plot(1:1:size(xi,2), xi(1:3,:)', 'LineWidth', 1.1)
    grid on
    xlabel('Segment Nr. $a$', 'Interpreter', 'latex')
    ylabel('Discrete deformation gradient $\xi_a$ (rot.)', 'Interpreter', 'latex')
    legend({'$\xi_{a,1}$', '$\xi_{a,2}$', '$\xi_{a,3}$'}, 'interpreter', 'latex');
    title('Discrete Deformation Gradients (Rot.)', 'Interpreter', 'latex');

    ax.LineStyleOrder = ["-o", "--+", ":x"];
    ax.LineStyleCyclingMethod = "withcolor";
    xlim([1, nSeg]);


    % Plot discrete deformations (translational)
    ax = nexttile;

    plot(1:1:size(xi,2), xi(4:6,:)', 'LineWidth', 1.1)
    grid on
    xlabel('Segment Nr. $a$', 'Interpreter', 'latex')
    ylabel('Discrete deformation gradient $\xi_a$ (transl.)', 'Interpreter', 'latex')
    legend({'$\xi_{a,4}$', '$\xi_{a,5}$', '$\xi_{a,6}$'}, 'interpreter', 'latex');
    title('Discrete Deformation Gradients (Transl.)', 'Interpreter', 'latex');
    
    ax.LineStyleOrder = ["-o", "--+", ":x"];
    ax.LineStyleCyclingMethod = "withcolor";
    xlim([1, nSeg]);

end
