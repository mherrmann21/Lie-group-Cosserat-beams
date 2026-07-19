function pars = beamParams_LLA11_steelString(options)
    %% Beam parameters for the "steel string" example from [LLA11]
    % with circular cross-section
    %
    %
    % Maximilian Herrmann
    % Chair of Automatic Control
    % TUM School of Engineering and Design
    % Technical University of Munich

    arguments
        % Individual Value for the beam length (m)
        options.L (1,1) double = 1;

        % Individual Value for Young's Modulus (N/m^2)
        % Default Value: Steel
        % https://de.wikipedia.org/wiki/Elastizit%C3%A4tsmodul#Typische_Zahlenwerte
        options.E (1,1) double = 2.1e11;

        % Cross-section radius, see [LLA11, p.307]
        options.radius (1,1) double = 1e-3;
    end

    % Class Instance
    pars = beamParams;


    %% Beam Parameters

    % Beam length (m)
    % Default value set above in arguments block
    pars.L = options.L;


    %%% Beam Geometry
    % with circular cross-section

    radius = options.radius;

    % Cross-section geometry
    % H/W corresponds to the diameter of the circular cross-section
    pars.geom.H = 2*radius;
    pars.geom.W = 2*radius;
    pars.geom.A = radius^2 * pi;

    % Compute second moments of inertia (about x and y axes of the body-fixed
    % coordinate systems)
    % https://en.wikipedia.org/wiki/List_of_second_moments_of_area
    pars.geom.I_x = pi/4 * radius^4;
    pars.geom.I_y = pi/4 * radius^4;

    % Polar moment of inertia
    pars.geom.J_P = pi/2 * radius^4;


    %%% Beam Material
    % Material: Steel

    % Density (kg/m^3)
    pars.mat.rho = 7.85e3;

    % Young's modulus (N/m^2)
    % https://de.wikipedia.org/wiki/Elastizit%C3%A4tsmodul#Typische_Zahlenwerte
    % Default value set above in arguments block
    pars.mat.E = options.E;

    % Poisson's ratio
    % https://en.wikipedia.org/wiki/Poisson%27s_ratio#Poisson's_ratio_values_for_different_materials
    pars.mat.nu = 0.2;

    %%% Dissipation coefficients
    pars.d = 0;

    pars = pars.computeParams;


    %% Display some key properties

    fprintf('Beam length-to-height ratio L/H: %d\n', pars.L/pars.geom.H);
end
