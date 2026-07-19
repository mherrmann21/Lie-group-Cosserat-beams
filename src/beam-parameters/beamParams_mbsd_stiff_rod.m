function pars = beamParams_mbsd_stiff_rod()
    %% Beam parameters for the mbsd paper
    %
    % Maximilian Herrmann
    % Chair of Automatic Control
    % TUM School of Engineering and Design
    % Technical University of Munich

    % Class Instance
    pars = beamParams;


    %% Beam Parameters

    % Beam length (m)
    pars.L = 1;

    %%% Beam Geometry
    % with circular cross-section

    % Cross-section radius
    radius = 2e-3;

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
    % Material: high-strength aluminum

    % Density (kg/m^3)
    pars.mat.rho = 2.85e3;

    % Young's modulus (N/m^2)
    % https://de.wikipedia.org/wiki/Elastizit%C3%A4tsmodul#Typische_Zahlenwerte
    % Default value set above in arguments block
    pars.mat.E = 7.2E10;

    % Poisson's ratio
    % https://en.wikipedia.org/wiki/Poisson%27s_ratio#Poisson's_ratio_values_for_different_materials
    pars.mat.nu = 0.35;

    %%% Dissipation coefficients
    pars.d = 0;

    pars = pars.computeParams;


    %% Display some key properties

    fprintf('Beam length-to-height ratio L/H: %d\n', pars.L/pars.geom.H);
end
