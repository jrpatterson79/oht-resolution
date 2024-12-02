%% Cleanup

clear; close all; clc;

%% Specify Directory
addpath(genpath('/.../.../'))

%% Describe the setup for the forward models
%PARAMETERS: Visualization and saving options
case_name = '2_freq';

%PARAMETERS: Associated with domain and testing setup
%Specify domain
% [domain] = equigrid_setup(x_disc,y_disc);
xmin = -100; xmax = -xmin; dx = 2;
ymin = xmin; ymax = xmax; dy = dx;

domain = struct('x',[],'y',[],'z',[]);
domain.x = [xmin : dx : xmax];
domain.y = [ymin : dy : ymax];
domain.z = [0 1];

%Specify boundary types and boundary values (x / y constant head
%boundaries, no flux in z)
bdry_types = [1; 1; 1; 1; 0; 0];
bdry_vals = zeros(6,1);
bdry_L = 1e-5*ones(6,1);
bdrys = struct('types',bdry_types,'vals',bdry_vals,'leaks',bdry_L);

% Ideal 9-spot Well Geometry
well_locs = [-20 -20; ...
             -20 0; ...
             -20 20; ...
             0 -20; ...
             0 0; ...
             0 20; ...
             20 -20; ...
             20 0; ...
             20 20 ...
            ];

num_wells = size(well_locs,1);

% List defining each observation. Columns are:
% pumping angular frequency, pumping well, Q_max, and observation well
% NOTE: The model will run the fastest if this list is sorted by angular
% frequency and pumping well (the first two columns)
%
V = 0.01;
% Q_max = 7e-5;

% P = 100;
P = [10 1600];
% P = [10 50 300 1600];
% P = [10 50 90 300 700 1600];
% P = [10 30 50 90 200 300 700 1600];

test_list = [];
tseries_length = [];
for i = 1:1:numel(P)
    for j = 1:1: num_wells
        for k = j+1:1:num_wells
            test_list = [test_list; ...
                         (2*pi)/P(i) j (V*pi)/P(i) k ...
                         ];
            tseries_length = [tseries_length;...
                5*P(i)];
        end
    end
end

% PARAMETERS: Associated with creating test case for fields

% Geostatistical parameters
lnK_mean = -9.2;   % Mean log(hydraulic conductivity (m/s))
lnSs_mean = -11.2; % Mean log(specific storage (1/m))

lnK_var = 2;    % Variance log(hydraulic conductivity (m/s))
lnSs_var = 0.4; % Variance log(specific storage (1/m))

corr_x = 20; % x correlation length (m)
corr_y = 20; % y correlation length (m)

% PARAMETERS: Associated with inversion setup
lnK_guess = -9;    % Initial guess hydraulic conductivity
lnSs_guess = -11;  % Initial guess specific storage

lnK_var_est = 1;    % Estimated variance - hydraulic conductivity
lnSs_var_est = 0.1; % Estimated variance - specific storage

%% Create the "official" input files needed by the steady periodic model.

%This step should automatically create the input files needed to run all of
%the steady periodic models.
[experiment] = OHT_create_inputs(well_locs,test_list,domain);

%% Calculate counts (used often for plotting and other functionality)
% Calculates number of observations, size of grid, and sets up a
% meshgridded set of points to do plotting of results.

num_omegas = size(experiment,1);
num_obs = size(test_list,1);
num_x = numel(domain.x) - 1;
num_y = numel(domain.y) - 1;
num_cells = num_x*num_y;

%% Setup grid of cell centers for plotting and geostatistical setups
[coords, cgrid] = plaid_cellcenter_coord(domain);

%% For synthetic problem - true parameter field statistics
num_relz = 100;
parfor a = 1 : num_relz
    randn('state', a)
    %Geostatistical Case (Exponential variogram)
    %TODO: Allow choice of different variogram models?
    distmat_row = dimdist(coords(1,:),coords);
    corr_row = exp(-(...
        (distmat_row(:,:,1)./corr_x).^2 + ...
        (distmat_row(:,:,2)./corr_y).^2).^.5);
    [corr_relz] = toepmat_vector_math(corr_row,'r',[],2,[num_y num_x]);
    lnK_true = corr_relz(:,1).*sqrt(lnK_var) + lnK_mean;
    lnSs_true = corr_relz(:,1) .* sqrt(lnSs_var) + lnSs_mean;
    params_true = [lnK_true; lnSs_true];

    lnK_range = [min(lnK_true) max(lnK_true)];
    lnSs_range = [min(lnSs_true) max(lnSs_true)];

    %% Perform all model runs to generate data

    % Steady_periodic_multimodel_run performs the model runs to generate all
    % observations.
    % IMPORTANT:
    % The first output is simulated observations, output as
    % follows: [A_obs(1); A_obs(2); A_obs(3) ... A_obs(numobs); B_obs(1) ...
    % B_obs(numobs)]
    %
    % NOTE: B in this case is the imaginary part of the phasor. Note that this
    % is not the same as the coefficient in front of the sine (which would be
    % negative B)
    %
    % The second output is the full phasor field for each pumping test.
    % TODO: May change this so that it is the given pumping test for each
    % observation. Also, need to add description of structure.
    %
    % H is the sensitivity matrix. Each row represents a different observation,
    % organized the same way as sim_obs. Each column represents a different
    % parameter sensitivity - the first set of columns is with respect to each
    % K value, and the second set of columns is with respect to each Ss value.
    % K and Ss grid cells are stepped through in sequential order by y, then x,
    % then z, as would be done by meshgrid.
    h_function = @(params) OHT_run_distribKSs(params, ...
        domain,bdrys,experiment,1);

    Phi_function = @(params) OHT_run_distribKSs(params, ...
        domain,bdrys,experiment,2);

    H_tilde_func = @(params) OHT_run_distribKSs(params, ...
        domain,bdrys,experiment,3);

    %Generate example simulated data
    tic
    sim_phasor = h_function(params_true);
    toc

    %Generate phasor field variable
    tic
    Phi_true = Phi_function(params_true);
    toc

    params_init = [lnK_guess*ones(num_cells,1); lnSs_guess*ones(num_cells,1)];

    %Generate example sensitivity map given initial homogeneous parameter
    %guesses
    tic
    H_adj = H_tilde_func(params_init);
    toc

    %% Add noise to synthetic data timeseries
    %Parameters used in creating synthetic time-series data and associated
    %noise
    tseries_sd = 0.5e-3;  % Head time-series error standard deviation (m)
    tseries_step = 1/125; % Time-series time step (125 Hz sampling frequency)

    obs_phasor_vec = zeros(2*num_obs,1);
    obs_phasor_covar = zeros(2*num_obs,2*num_obs);
    for o = 1:1:num_obs

        t = [0:tseries_step:tseries_length(o)]';
        num_t = numel(t);
        %Creation of noisy time series data
        om = test_list(o,1);
        cos_ts = cos(om*t);
        sin_ts = sin(om*t);
        M = [cos_ts -sin_ts];
        tseries_noisy = sim_phasor(o)*cos_ts - sim_phasor(o+num_obs)*sin_ts + randn(num_t,1)*tseries_sd;

        %Extraction of phasor and error estimates
        obs_phasor = (M'*M)\(M'*tseries_noisy);
        tseries_err_est = (sum((tseries_noisy - M*obs_phasor).^2)./(num_t - 2))^.5;
        disp(['Noise-free phasor: ', num2str(sim_phasor(o)), ', ', num2str(sim_phasor(o+num_obs))])
        disp(['Estimated phasor: ', num2str(obs_phasor(1)), ', ', num2str(obs_phasor(2))])
        disp(['Time Series Noise Estimate: ', num2str(tseries_err_est)])
        phasor_covar_est = (tseries_err_est.^2*inv(M'*M));

        %Placements of phasors into observation vector, and covariance values
        %(calculated via linear covariance propagation) into covariance matrix.
        obs_phasor_vec(o) = obs_phasor(1);
        obs_phasor_vec(o+num_obs) = obs_phasor(2);
        obs_phasor_covar(o,o) = phasor_covar_est(1,1);
        obs_phasor_covar(o,o+num_obs) = phasor_covar_est(1,2);
        obs_phasor_covar(o+num_obs,o) = phasor_covar_est(2,1);
        obs_phasor_covar(o+num_obs,o+num_obs) = phasor_covar_est(2,2);
    end

    %% Prior setup
    %TODO: Allow choice of different variogram models?
    % Assumed linear variogram prior parameter covariance
    distmat_row = dimdist(coords(1,:),coords);
    max_dist = max(max((distmat_row(:,:,1).^2 + distmat_row(:,:,2).^2).^.5));
    slope = 1/max_dist;

    Q_row_corr = slope .*...
                (max_dist - ((distmat_row(:,:,1).^2 + ...
                 distmat_row(:,:,2).^2).^.5));

    %% Inversion
    R = obs_phasor_covar; % Data error covariance matrix
    beta_init = [lnK_guess; lnSs_guess]; % Initial parameter guess
    X = [ones(num_cells,1) zeros(num_cells,1); zeros(num_cells,1) ones(num_cells,1)]; % Parameter drift matrix
    y = obs_phasor_vec;
    coords_range = (coords(:,1) > -40 & coords(:,1) < 40 & coords(:,2) > -40 & coords(:,2) < 40);

    Q_row_K = Q_row_corr .* lnK_var_est;
    Q_row_Ss = Q_row_corr .* lnSs_var_est;
    Qprod_func = @(invec) covar_product_K_Ss(Q_row_K, Q_row_Ss,invec,num_x,num_y);

    tic;
    [params_best, beta_best, H_local, NLAP] = ql_geostat_LM_inv(y,params_init,beta_init,X,R,Qprod_func,h_function,H_tilde_func);
    inv_time = toc;

    negloglike_func = @(s,beta) negloglike_eval(y,X,s,beta,Q,R,h_function);
    y_obs_mat(:,a) = y;
    params_true_mat(:,a) = params_true;

    params_best_mat(:,a) = params_best; % Optimal parameters
    beta_best_vec(a,:) = beta_best; % Optimal mean parameters
end

%% Save output
save_dir = '/.../.../'; % Directory to save output .mat file
save_name = ['stoch_' case_name '.mat'];
save([save_dir save_name])