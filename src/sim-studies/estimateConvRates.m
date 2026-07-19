function [convTime, convSpace] = estimateConvRates(values, E, h, nSeg, modelNames)
    arguments
        values      (:,:,:,:) double
        E           (:,1)
        h           (:,1) double
        nSeg        (:,1) double
        modelNames  (:,1) string
    end

    % Set up fittype and options
    ft = fittype( 'power1' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';

    % Levenberg-Marquardt and tighter tolerances seem to give much better
    % fit
    opts.Algorithm = 'Levenberg-Marquardt';
    opts.TolFun = 1e-10;
    opts.TolX   = 1e-10;
    opts.StartPoint = [100, 1.2];

    % Convergence over time
    convTime.rsquare = nan(length(E), length(nSeg), length(modelNames));
    convTime.a = nan(length(E), length(nSeg), length(modelNames));
    convTime.b = nan(length(E), length(nSeg), length(modelNames));

    for iMat = 1:length(E)
        for iSeg = 1:length(nSeg)
            for iSim = 1:length(modelNames)

                % Check if we have enough data for the fit
                % (at least 3 non-nan values)
                if nnz(~isnan(values(iMat, :, iSeg, iSim) )) > 2

                    % Prepare fit (code generated from matlab curve fitter app)
                    % Note: Use the quadratic time step as weights to get
                    % better fit due to the exponentially small values of h
                    [xData, yData, weights] = prepareCurveData( ...
                        h, values(iMat, :, iSeg, iSim)', (1./h).^2 ...
                        );

                    opts.Weights = weights;

                    % Fit model to data
                    [fitresult, gof] = fit( xData, yData, ft, opts );

                    % Get error and coefficients
                    convTime.rsquare(iMat, iSeg, iSim) = gof.rsquare;
                    convTime.a(iMat, iSeg, iSim) = fitresult.a;
                    convTime.b(iMat, iSeg, iSim) = fitresult.b;
                end
            end
        end
    end


    % Convergence over space
    convSpace.rsquare = nan(length(E), length(h), length(modelNames));
    convSpace.a = nan(length(E), length(h), length(modelNames));
    convSpace.b = nan(length(E), length(h), length(modelNames));

    for iMat = 1:length(E)
        for ih = 1:length(h)
            for iSim = 1:length(modelNames)

                % Check if we have enough data for the fit
                % (at least 3 non-nan values)
                if nnz(~isnan(values(iMat, ih, :, iSim) )) > 2

                    % Prepare fit (code generated from matlab curve fitter app)
                    [xData, yData, weights] = prepareCurveData( ...
                        1./nSeg, ...
                        squeeze( values(iMat, ih, :, iSim) ), ...
                        nSeg.^2 ...
                        );

                    opts.Weights = weights;

                    % Fit model to data
                    [fitresult, gof] = fit( xData, yData, ft, opts );

                    % Get error and coefficients
                    convSpace.rsquare(iMat, ih, iSim) = gof.rsquare;
                    convSpace.a(iMat, ih, iSim) = fitresult.a;
                    convSpace.b(iMat, ih, iSim) = fitresult.b;
                end
            end
        end
    end

end