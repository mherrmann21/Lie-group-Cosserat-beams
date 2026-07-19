function animateBeamSimRes(simData, params, SAVE_MOVIE, fileName, opts)
    %% Animate the simulation results of a beam simulation
    % Todo:
    % * Proper description and input definition
    % * Proper definition of the simRes struct and its fields
    % * Option to specify video file name

    arguments
        simData     (1,1) beamSimData
        params      (1,1) beamParams

        % Save movie to file?
        SAVE_MOVIE  (1,1) logical

        % File name/path for the saved movie
        fileName    (1,1) string

        opts.frameRate (1,1) double = 30;
    end

    % Disable warning about video autput size
    warning('off', 'MATLAB:audiovideo:VideoWriter:mp4FramePadded');

    tout = simData.tout;
    %xout = simRes.xout;

    nNodes = size(simData.R, 3);
    nSeg = nNodes-1;

    beamW = params.geom.W;
    beamH = params.geom.H;


    %% Interpolate results at fixed sampling rate and spatial resolution

    tSample = 1/opts.frameRate; % For 30FPS video
    tQuery = (tout(1):tSample:tout(end))';
    gInput = SE3Matrix(simData.R,simData.x);

    % Use linear interpolation in time
    gOutT = interpn( ...
        1:4, 1:4, 1:size(gInput,3), simData.tout, ...
        gInput, ...
        1:4, 1:4, 1:size(gInput,3), tQuery);

    xiInput = interpn( ...
        1:6, 1:size(gInput,3)-1, simData.tout, ...
        simData.xi, ...
        1:6, 1:size(gInput,3)-1, tQuery ...
        );

    % Spatial interpolation on SE3

    % Spatial query points (normalize beam length to 1)
    sInput = 0:(params.L/nSeg):params.L;
    lInterp = 0.01;
    sQuery = 0:lInterp:params.L;
    gQuery = interpSimResSpaceSE3(gOutT,xiInput,sInput,sQuery);


    %% Animate Beam

    % Disable warning for badly conditioned SE3 matrices (due to linear
    % interpolation)
    warning('off', 'MATLAB:hg:DiceyTransformMatrix');

    [fig,ax] = init3Dplot( 'Name', 'SimAnim' );

    % Make background color white for nice videos
    fig.Color = [1,1,1];

    ax.TickLabelInterpreter = 'latex';

    if ~strcmp(fig.WindowStyle, 'docked')
        fig.WindowState = 'maximized';
    end

    xlabel('$x$ / m', 'Interpreter', 'latex'); 
    ylabel('$y$ / m', 'Interpreter', 'latex');  
    zlabel('$z$ / m', 'Interpreter', 'latex'); 
    axis equal

    beamVis = elasticBeam(gQuery(:,:,:,1), ...
        'showFrames', false,...
        'height', beamH, 'width', beamW ...
        );

    % Add text for current time step
    % Place reasonably high above the beam
    zPos = 0.5*beamH + 0.05;
    textTime = text(0,0,0.3,'', 'Interpreter', 'latex');

    %title('Simulation Results', 'Interpreter', 'latex');

    % Extra margin
    sc = beamH+0.05;

    xlim([ ...
        min(simData.x(1,:,:),[], 'all')-sc, ...
        max(simData.x(1,:,:),[], 'all')+sc] ...
        );
    ylim([ ...
        min(simData.x(2,:,:),[], 'all')-sc, ...
        max(simData.x(2,:,:),[], 'all')+sc] ...
        );
    zlim([ ...
        min(simData.x(3,:,:),[], 'all')-sc, ...
        max(max(simData.x(3,:,:),[], 'all')+sc, zPos)] ... % Consider position of the time label
        );
    drawnow;

    % Preallocate animation frames struct 
    % (taken from getframe documentation)
    animFrame(length(1:length(tQuery))) = struct('cdata',[],'colormap',[]);

    for iStep = 1:length(tQuery)

        % Check if figure is closed
        if ~isvalid(fig)
            return
        else
            textTime.String = sprintf('$t$ = %.3fs', tQuery(iStep));

            beamVis.updateConfiguration(gQuery(:,:,:,iStep));
            drawnow;

            if isvalid(fig) && SAVE_MOVIE
                try
                    animFrame(iStep) = getframe(fig);
                catch
                    return;
                end
            end
        end


    end

    %% Write to video
    % (code from MATLAB docs)
    if SAVE_MOVIE
        disp('Saving as Video...')

        if isvalid(fig)
            v = VideoWriter(fileName, 'MPEG-4');
            v.Quality = 100;
            open(v);
            for iFrame = 1:length(animFrame)
                writeVideo(v,animFrame(iFrame));
            end

            close(v);
        end
    end
end
