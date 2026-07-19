function pars = beamParams_outOfPlaneLoad()
    %% Sample parameters for beam simulation
    % Parameters for the "out of plane load test" from literature,
    % e.g., [How+18]
    %
    % Maximilian Herrmann
    % Chair of Automatic Control
    % TUM School of Engineering and Design
    % Technical University of Munich

    % Class instance
    pars = beamParams;

    % Beam length (m)
    % Compute length from desired beam shape
    radius = 100;
    pars.L = 2*radius*pi/8;


    %%% Beam Geometry
    % with rectangular cross-section

    % Cross-section geometry
    pars.geom.H = 1;
    pars.geom.W = 1;
    pars.geom.A = pars.geom.H * pars.geom.W;

    % Compute second moments of inertia (about x and y axes of the body-fixed
    % coordinate systems)
    % https://en.wikipedia.org/wiki/List_of_second_moments_of_area
    pars.geom.I_x =  pars.geom.H^4 / 12;
    pars.geom.I_y =  pars.geom.W^4 / 12;
    pars.geom.J_P = pars.geom.I_x + pars.geom.I_y;


    %%% Beam Material
    pars.mat.E = 1e7;
    pars.mat.nu = 0;

    % From the LLA12 rubber-rod example (not given and not relevant for the
    % static simulation)
    pars.mat.rho = 1e0;

    %%% Dissipation coefficients
    pars.d = 0;

    pars = pars.computeParams;

end
