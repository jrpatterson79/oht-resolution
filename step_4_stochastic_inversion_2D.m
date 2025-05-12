% Oscillatory Hydraulic Tomography (OHT) Linear Resolution Analysis

% This code conducts stochastic geostatistical inverse modeling. In a parallelized loop, this code generates 100 heterogeneous transmissivity and storativity realizations and then conducts inverse modeling for single- and multi-frequency OHT. The code outputs .mat files with inversion results that are used for plotting in step_5_plot_stochastic_results_2D.m.

% The user chooses single- or multi-frequency OHT testing by commenting / uncommenting lines 64-70.

% Unmodified, this code will reproduce the analysis presented in:
% Patterson, J. R., & Cardiff, M. (2025). Multi‐frequency oscillatory hydraulic tomography improves heterogeneity imaging and resolution and reduces uncertainty. Water Resources Research, 61, e2024WR039606. https://doi.org/10.1029/2024WR039606

% This code requires OHT3DINV v.0.16.0 which can be downloaded at https://github.com/wischydro-cardiff/oscillatory-tomography.git

% Code developed by Jeremy Patterson
% Created: March 2022; Updated May 2025

%% Clean Environment
clear; close all; clc;

%% Specify Directory for OHT3DINV
addpath(genpath('/.../.../.../'))

%% Describe the setup for the forward models
%PARAMETERS: Visualization and saving options
case_name = '4_freq';

%PARAMETERS: Associated with domain and testing setup
%Specify domain
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

% P = 10;
% P = 100;
% P = 1600;
% P = [10 1600];
P = [10 50 300 1600];
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
lnT_mean = -9.2;   % Mean log(transmissivity (m^2/s))
lnS_mean = -11.2; % Mean log(storativity (-))

lnT_var = 2;    % Variance log(transmissivity (m^2/s))
lnS_var = 0.5; % Variance log(storativity (-))

corr_x = 20; % x correlation length (m)
corr_y = 20; % y correlation length (m)


% PARAMETERS: Associated with inversion setup
lnT_guess = -9;    % Initial guess transmissivity
lnS_guess = -11;  % Initial guess storativity

lnT_var_est = 1;    % Estimated variance - transmissivity
lnS_var_est = 0.5; % Estimated variance - storativity

%% Create the "official" input files needed by the steady periodic model.

%This step automatically creates the input files needed to run all of
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
    % Generate heterogeneous parameter fields (Exponential variogram model)
    distmat_row = dimdist(coords(1,:),coords);
    corr_row = exp(-(...
        (distmat_row(:,:,1)./corr_x).^2 + ...
        (distmat_row(:,:,2)./corr_y).^2).^.5);
    [corr_relz] = toepmat_vector_math(corr_row,'r',[],2,[num_y num_x]);
    lnT_true = corr_relz(:,1).*sqrt(lnT_var) + lnT_mean;
    lnS_true = corr_relz(:,1) .* sqrt(lnS_var) + lnS_mean;
    params_true = [lnT_true; lnS_true];

    lnT_range = [min(lnT_true) max(lnT_true)];
    lnS_range = [min(lnS_true) max(lnS_true)];

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

    params_init = [lnT_guess*ones(num_cells,1); lnS_guess*ones(num_cells,1)];

    %Generate example sensitivity map given initial homogeneous parameter
    %guesses
    tic
    H_adj = H_tilde_func(params_init);
    toc

    %% Add noise to synthetic data timeseries
    %Parameters used in creating synthetic time-series data and associated
    %noise
    tseries_sd = 0.5e-3;
    tseries_step = 1/125;

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

    %% Prior setup - Linear Variogram Model
    max_dist = max(max((distmat_row(:,:,1).^2 + distmat_row(:,:,2).^2).^.5));
    slope = 1/max_dist;

    Q_row_corr = slope .*...
                (max_dist - ((distmat_row(:,:,1).^2 + ...
                 distmat_row(:,:,2).^2).^.5));

    %% Inversion
    R = obs_phasor_covar;
    beta_init = [lnT_guess; lnS_guess];
    X = [ones(num_cells,1) zeros(num_cells,1); zeros(num_cells,1) ones(num_cells,1)];
    y = obs_phasor_vec;
    coords_range = (coords(:,1) > -40 & coords(:,1) < 40 & coords(:,2) > -40 & coords(:,2) < 40);

    Q_row_K = Q_row_corr .* lnT_var_est;
    Q_row_Ss = Q_row_corr .* lnS_var_est;
    Qprod_func = @(invec) covar_product_K_Ss(Q_row_K, Q_row_Ss,invec,num_x,num_y);

    tic;
    [params_best, beta_best, H_local, NLAP] = ql_geostat_LM_inv(y,params_init,beta_init,X,R,Qprod_func,h_function,H_tilde_func);
    inv_time = toc;

    negloglike_func = @(s,beta) negloglike_eval(y,X,s,beta,Q,R,h_function);
    y_obs_mat(:,a) = y;
    params_true_mat(:,a) = params_true;

    params_best_mat(:,a) = params_best;
    beta_best_vec(a,:) = beta_best;
end

%% Save output
save_dir = '/.../.../.../';
save_name = ['stoch_' case_name '.mat'];
save([save_dir save_name])