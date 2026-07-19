function figHandles = plotQuantity2d(values, dimInfo, figInfo, convTime, convSpace)
    arguments
        values      (:,:,:,:)   double
        dimInfo     (:,1)       struct
        figInfo     (1,1)       struct

        % Struct with convergence rate info over time; can be empty
        % Dimensions: (E,nSeg,integrators)
        convTime    (:,:)       struct

        % Struct with convergence rate info over space; can be empty
        % Dimensions: (E,h,integrators)
        convSpace   (:,:)       struct
    end

    % Expected Data dimensions: (E, h, nSeg, integrators)

    %% 1. Value over time step, with all integrators in one plot
    % X-Axis  (Data dim 1): h
    % Lines   (Data dim 2): Integrators
    % Rows    (Data dim 3): E
    % Columns (Data dim 4): nSeg

    dimOrder = [2,4,1,3];

    if ~isempty(convTime)
        % convLines.a = mean(convTime.a, 3, "omitnan" );
        % convLines.b = mean(convTime.b, 3, "omitnan" );
        convLines.a = permute(convTime.a, [1,2,3]);
        convLines.b = permute(convTime.b, [1,2,3]);
    else
        convLines = [];
    end

    [figHandles(1), ~] = plot4DValueGrid( ...
        permute(values, dimOrder),  dimInfo(dimOrder),  figInfo, ...
        convLines ...
        );


    %% 2. Value over time step, with all segment numbers in one plot
    % X-Axis  (Data dim 1): Time step
    % Lines   (Data dim 2): nSeg
    % Rows    (Data dim 3): E
    % Columns (Data dim 4): Integrators

    dimOrder = [2,3,1,4];

    if ~isempty(convTime)
        % convLines.a = squeeze(mean(convTime.a, 2, "omitnan" ));
        % convLines.b = squeeze(mean(convTime.b, 2, "omitnan" ));
        convLines.a = permute(convTime.a, [1,3,2]);
        convLines.b = permute(convTime.b, [1,3,2]);
    else
        convLines = [];
    end

    [figHandles(2), ~] = plot4DValueGrid( ...
        permute(values, dimOrder),  dimInfo(dimOrder),  figInfo, ...
        convLines ...
        );


    %% 3. Value over Nr. of Segments, with all integrators in one plot
    % X-Axis  (Data dim 1): nSeg
    % Lines   (Data dim 2): Integrators
    % Rows    (Data dim 3): E
    % Columns (Data dim 4): h

    dimOrder = [3,4,1,2];

    if ~isempty(convSpace)
        % convLines.a = squeeze(mean(convSpace.a, 3, "omitnan" ));
        % convLines.b = squeeze(mean(convSpace.b, 3, "omitnan" ));
        convLines.a = permute(convSpace.a, [1,2,3]);
        convLines.b = permute(convSpace.b, [1,2,3]);
    else
        convLines = [];
    end

    [figHandles(3), ~] = plot4DValueGrid( ...
        permute(values, dimOrder),  dimInfo(dimOrder),  figInfo, ...
        convLines ...
        );


    %% 4. Value over Nr. of Segments, with all time steps in one plot
    % X-Axis  (Data dim 1): nSeg
    % Lines   (Data dim 2): h
    % Rows    (Data dim 3): E
    % Columns (Data dim 4): Integrators

    dimOrder = [3,2,1,4];

    if ~isempty(convSpace)
        % convLines.a = squeeze(mean(convSpace.a, 2, "omitnan" ));
        % convLines.b = squeeze(mean(convSpace.b, 2, "omitnan" ));
        convLines.a = permute(convSpace.a, [1,3,2]);
        convLines.b = permute(convSpace.b, [1,3,2]);
    else
        convLines = [];
    end

    [figHandles(4), ~] = plot4DValueGrid( ...
        permute(values, dimOrder),  dimInfo(dimOrder),  figInfo, ...
        convLines ...
        );


    %% 5. Value over Young's Modulus
    % X-Axis  (Data dim 1): E
    % Lines   (Data dim 2): Integrators
    % Rows    (Data dim 3): nSeg
    % Columns (Data dim 4): h

    if length(dimInfo(1).xValues) > 1

        dimOrder = [1,4,3,2];

        [figHandles(5), ~] = plot4DValueGrid( ...
            permute(values, dimOrder),  dimInfo(dimOrder),  figInfo, [] ...
            );

    end

end