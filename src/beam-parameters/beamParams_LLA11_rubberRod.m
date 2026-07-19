function pars = beamParams_LLA11_rubberRod(options)
    %% Beam parameters for the "rubber rod" example from [LLA11]
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
        % See [LLA11, p.301]
        options.E (1,1) double = 5e6;
    end

    % Class Instance
    pars = beamParams;


    %% Beam Parameters

    % Beam length (m)
    % Default value set above in arguments block
    pars.L = options.L;


    %%% Beam Geometry
    % with circular cross-section

    % Cross-section radius, see [LLA11, p.301]
    radius = 5e-3;

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
    % Material: Rubber, see [LLA11, p.301]

    % Density (kg/m^3)
    pars.mat.rho = 1.1e3;

    % Young's modulus (N/m^2)
    % Default value set above in arguments block
    pars.mat.E = options.E;

    % Poisson's ratio
    % Actually 0.5 in the paper! 0.49 to be able to compute CFL limit
    %pars.mat.nu = 0.49;
    pars.mat.nu = 0.5;


    %%% Dissipation coefficients
    pars.d = 0;

    pars = pars.computeParams;
    

    %% Display some key properties

    fprintf('Beam length-to-height ratio L/H: %d\n', pars.L/pars.geom.H);
end
