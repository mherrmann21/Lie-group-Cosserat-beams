function hCFL = checkCFLLimit(beamPars, nSeg, options)
    % Check the beam CFL limit for a given discretization and material
    % Formulas after [Dem+14], Remark 5.2, p.30

    arguments
        beamPars
        nSeg
        options.h = nan;
    end

    %% CFL Limit

    nu = beamPars.mat.nu;
    E = beamPars.mat.E;
    G = beamPars.mat.G;

    % Lamé constants
    % https://de.wikipedia.org/wiki/Lam%C3%A9-Konstanten
    % lambda = nu/(1-2*nu) * 1 / (1+nu)*E
    lambda = G*(E-2*G)  / ( 3*G - E );
    mu = beamPars.mat.G;

    d = beamPars.L / nSeg;
    c = sqrt( (lambda +2*mu) / beamPars.mat.rho  );

    hCFL = d / (10*c);

    fprintf('   Maximum CFL time step: %e = 2^%.2f s\n', hCFL, log2(hCFL));

    if ~isnan(options.h)
        fprintf('   Used time step:        %e = 2^%.2f s\n', options.h, log2(options.h));
        if options.h > hCFL
            warning('Chosen time step (%e = 2^%.2f s) is above maximum timestep according to CFL limit (%e = 2^%.2f s)!\n', options.h, log2(options.h), hCFL, log2(hCFL) );
        end
    end

end
