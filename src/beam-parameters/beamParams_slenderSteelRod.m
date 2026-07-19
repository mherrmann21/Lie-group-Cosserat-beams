function pars = beamParams_slenderSteelRod(options)
    %% Beam parameters for a slender steel beam
    % with quadratic/rectangular cross-section
    %
    %
    % Maximilian Herrmann
    % Chair of Automatic Control
    % TUM School of Engineering and Design
    % Technical University of Munich
    %
    % The beam length and Young's modulus can be overridden with name-value
    % arguments.

    arguments
        % Individual Value for the beam length (m)
        options.L (1,1) double = 1;

        % Individual Value for Young's Modulus (N/m^2)
        % Default Value: Steel
        % https://de.wikipedia.org/wiki/Elastizit%C3%A4tsmodul#Typische_Zahlenwerte
        options.E (1,1) double = 2.1e11;
    end

    % Class Instance
    pars = beamParams;


    %% Beam Parameters

    % Beam length (m)
    % Default value set above in arguments block
    pars.L = options.L;


    %%% Beam Geometry
    % with rectangular cross-section

    % Cross-section geometry
    pars.geom.H = 0.01;
    pars.geom.W = 0.01;
    pars.geom.A = pars.geom.H * pars.geom.W;

    % Compute second moments of inertia (about x and y axes of the body-fixed
    % coordinate systems)
    % https://en.wikipedia.org/wiki/List_of_second_moments_of_area
    pars.geom.I_x = pars.geom.A * pars.geom.H^2 / 12;
    pars.geom.I_y = pars.geom.A * pars.geom.W^2 / 12;

    % Polar moment of inertia
    pars.geom.J_P = pars.geom.I_x + pars.geom.I_y;


    %%% Beam Material
    % Material: Steel

    % Density (kg/m^3)
    pars.mat.rho = 8050;

    % Young's modulus (N/m^2)
    % https://de.wikipedia.org/wiki/Elastizit%C3%A4tsmodul#Typische_Zahlenwerte
    % Default value set above in arguments block
    pars.mat.E = options.E;

    % Poisson's ratio
    % https://en.wikipedia.org/wiki/Poisson%27s_ratio#Poisson's_ratio_values_for_different_materials
    pars.mat.nu = 0.3;

    %%% Dissipation coefficients
    pars.d = 0;

    pars = pars.computeParams;

    
    %% Display some key properties

    fprintf('Beam length-to-height ratio L/H: %d\n', pars.L/pars.geom.H);
end
