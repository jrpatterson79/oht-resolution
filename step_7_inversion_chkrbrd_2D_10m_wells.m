% Oscillatory Hydraulic Tomography (OHT) Linear Resolution Analysis

% This code conducts checkerboard testing to assess resolution associated with single-frequency OHT with 10-m well separation. The code generates a checkerboard pattern of heterogeneous transmissivity and storativity then produces synthetic data. The generated data is used in geostatistical inverse modeling to determine how well the checkerboard is recovered. The code outputs .mat files with inversion results that are used for plotting in step_8_plot_chkrbrd_results_2D.m.

% The user chooses the pumping frequency on line 82. The user chooses the size of individual checkers on lines 106 and 107.

% Unmodified, this code will reproduce the analysis presented in:
% Patterson, J. R., & Cardiff, M. (2025). Multi‐frequency oscillatory hydraulic tomography improves heterogeneity imaging and resolution and reduces uncertainty. Water Resources Research, 61, e2024WR039606. https://doi.org/10.1029/2024WR039606

% This code requires OHT3DINV v.0.16.0 which can be downloaded at https://github.com/wischydro-cardiff/oscillatory-tomography

% Code developed by Jeremy Patterson
% Created: March 2022; Updated May 2025
%% Clean Environment
clear all; close all; clc;

%% Specify Directory to OHT3DINV
addpath(genpath('/.../.../.../'))

%% Describe the setup for the forward models

%PARAMETERS: Visualization and saving options
pause_length = 0.01;
case_name = '1_freq';

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
             -20 -10; ...
             -20 0; ...
             -20 10; ...
             -20 20; ...
             -10 -20; ...
             -10 -10; ...
             -10 0
             -10 10; ...
             -10 20; ...
             0 -20; ...
             0 -10; ...
             0 0; ...
             0 10; ...
             0 20; ...
             10 -20; ...
             10 -10; ...
             10 0
             10 10; ...
             10 20; ...
             20 -20; ...
             20 -10; ...
             20 0; ...
             20 10; ...
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

P = 1600;

test_list = [];
tseries_length = [];
for i = 1:1:numel(P)
    for j = 1:1: num_wells
        for k = j+1:1:num_wells
            test_list = [...
                test_list; ...
                (2*pi)/P(i) j (V*pi)/P(i) k ...
                ];
            tseries_length = [tseries_length;...
                              5*P(i)];
        end
    end
end

% PARAMETERS: Associated with creating test case for fields
%Used by both methods
lnT_mean = -9.2;   % Mean log(transmissivity (m^2/s))
lnS_mean = -11.2; % Mean log(storativity (-))

% Checkerboard geometry design
% Individual checker size
xl_check = 10; x_check_offset = 5;
yl_check = 10; y_check_offset = 2.5;
% Parameter deviation from mean value
lnT_jump = 0.5;
lnS_jump = 0.25;

% PARAMETERS: Associated with inversion setup
lnT_guess = -9;   % Transmissivity initial guess
lnS_guess = -11; % Storativity initial guess

lnT_var_est = 1;    % Variance log(transmissivity (m^2/s))
lnS_var_est = 0.5; % Variance log(storativity (-))

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

%% Visualize inputs

% Each time the angular frequency changes in "test_list", a new "omega
% group" is created. This is a set of models that get run together, which
% saves time.
%
% This block of code is just to make sure that all of the inputs for the
% steady periodic model have been created correctly. It will show:
% 1) In the command window, the current "omega group"
% 2) In Figure 1, the observation weighting for the current observation
% 3) In figure 2, the flux weighting for the given pumping test
%
% This will loop through all observations as you hit <Enter>, showing the
% model inputs for each individual model run that generates an observation.

fig_srcobs_weights = figure(1);
set(fig_srcobs_weights,'Position',[58 726 1100 600])
pause(1)
for i = 1:1:num_omegas
    disp(['Period = ', num2str(2*pi./experiment(i).omega)]);
    num_testobs = size(experiment(i).tests,1);
    for j = 1:1:num_testobs
        pump_loc = experiment(i).tests(j,1);
        obs_loc = experiment(i).tests(j,2);
        subplot(1,2,1)
        pumpmap = reshape(full(experiment(i).stims(:,pump_loc)),num_y,num_x);
        pm = pcolor(cgrid{1},cgrid{2},pumpmap);
        set(pm,'LineStyle','none')
        colorbar
        title(['Period group ', num2str(i), ', Pump weights, test/obs pair #', num2str(j)]);
        axis equal
        axis([xmin xmax ymin ymax])
        subplot(1,2,2)
        obsmap = reshape(full(experiment(i).obs(:,obs_loc)),num_y,num_x);
        om = pcolor(cgrid{1},cgrid{2},obsmap);
        set(om,'LineStyle','none')
        colorbar
        title(['Period group ', num2str(i), ', Observation weights, test/obs pair #', num2str(j)]);
        axis equal
        axis([xmin xmax ymin ymax])

        pause(pause_length)
    end
end

%% For synthetic problem - true parameter field statistics
% Checkerboard pattern
check_pattern = (sin(pi*(cgrid{1}-x_check_offset)/xl_check)).* ...
    (sin(pi*(cgrid{2}-y_check_offset)/yl_check));

lnT_true_grid = lnT_mean + check_pattern.*lnT_jump;
lnS_true_grid = lnS_mean + check_pattern.*lnS_jump;

lnT_true = reshape(lnT_true_grid,num_cells,1);
lnS_true = reshape(lnS_true_grid,num_cells,1);
params_true = ([lnT_true; lnS_true]);
    
lnT_range = [min(lnT_true) max(lnT_true)];
lnS_range = [min(lnS_true) max(lnS_true)];

figure
clf
subplot(1,2,1)
ax = gca;
p3 = pcolor(cgrid{1}, cgrid{2}, reshape(lnT_true, num_y, num_x));
p3.LineStyle = 'none';
hold on
plot(well_locs(:,1), well_locs(:,2), 'ko', 'LineWidth', 2,...
     'MarkerFaceColor', [192/255 16/255 0], 'MarkerSize', 9)
axis equal; axis square
axis([-50 50 -50 50])
ax.XTick = [-50:25:50]; 
ax.YTick = [-50:25:50];
xlabel('X direction (m)')
ylabel('Y direction (m)')
c = colorbar;
c.Label.String = 'ln(T [m^2/s])';
c.FontSize = 20;
caxis([lnT_range(1) lnT_range(2)])
ax.FontSize = 20;

subplot(1,2,2)
ax = gca;
p4 = pcolor(cgrid{1}, cgrid{2}, reshape(lnS_true, num_y, num_x));
p4.LineStyle = 'none';
hold on
plot(well_locs(:,1), well_locs(:,2), 'ko', 'LineWidth', 2,...
     'MarkerFaceColor', [192/255 16/255 0], 'MarkerSize', 9)
axis equal; axis square
axis([-50 50 -50 50])
ax.XTick = [-50:25:50]; 
ax.YTick = [-50:25:50];
xlabel('X direction (m)')
ylabel('Y direction (m)')
c = colorbar;
c.Label.String = 'ln(S [-])';
c.FontSize = 18;
caxis([lnS_range(1) lnS_range(2)])
ax.FontSize = 20;
set(gcf, 'Position', [100 100 2025 2025/2.75])
pause
close

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
    
    t = (0:tseries_step:tseries_length(o))';
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

    

%% Cross-checks - visualize results and sensitivities

%These cross-checks are used simply to make sure that the model appears to
%be running correctly. We calculate the amplitude and phase first (which
%are more intuitive than A and B), and then look at how these vary with
%space.
%
%The second set of plots looks at the sensitivity structure of the signal
%amplitude to K and Ss, which should look like "blobs" centered around the
%pumping and observation well for each observation.

A_obs = sim_phasor(1:num_obs);
B_obs = sim_phasor((num_obs+1):(2*num_obs));
amp = (A_obs.^2 + B_obs.^2).^.5;
phase = atan2(-B_obs,A_obs);

A_full = Phi_true(1:num_cells,:);
B_full = Phi_true((num_cells+1):(2*num_cells),:);
amp_full = (A_full.^2 + B_full.^2).^.5;
phase_full = atan2(-B_full,A_full);

H_AK = H_adj((1:num_obs),(1:num_cells));
H_BK = H_adj((num_obs+1):(2*num_obs),(1:num_cells));

H_ASs = H_adj((1:num_obs),((num_cells+1):(2*num_cells)));
H_BSs = H_adj((num_obs+1):(2*num_obs),((num_cells+1):(2*num_cells)));

% Check all phasor fields produced through numerical model running
fig_phaseamp = figure(3);
set(fig_phaseamp,'Position',[859 727 1100 600])
pause(1)

num_totalstims = size(amp_full,2);
for i = 1:1:num_totalstims
    amp_field = reshape(log(amp_full(:,i)),num_y,num_x);
    phase_field = reshape(phase_full(:,i),num_y,num_x);
    subplot(1,2,1);
    pc2 = pcolor(cgrid{1},cgrid{2},amp_field);
    set(pc2,'LineStyle','none')
    title(['Pumping test ', num2str(i), ', ln(Amplitude) field'])
    axis equal
    axis([xmin xmax ymin ymax])
    subplot(1,2,2);
    pc3 = pcolor(cgrid{1},cgrid{2},phase_field);
    set(pc3,'LineStyle','none')
    title(['Pumping test ', num2str(i), ', Phase field'])
    axis equal
    axis([xmin xmax ymin ymax])
    pause(pause_length)
end

for i = 1:1:num_obs
    ampK_sens = (A_obs(i)./amp(i)).*H_AK(i,:) + (B_obs(i)./amp(i)).*H_BK(i,:);
    ampK_sensmap = reshape(ampK_sens,num_y,num_x);
    figure
    clf
    ax = gca;
    pc4 = contourf(cgrid{1},cgrid{2}, ampK_sensmap, 10);
%     set(pc4,'LineStyle','none')
    hold on
    plot(well_locs(:,1), well_locs(:,2), 'ko', 'MarkerFaceColor', 'r', 'MarkerSize', 14)
    axis equal
    axis([-35 35 -25 45])
    xlabel('X Distance (m)')
    ylabel('Y Distance (m)')
    c = colorbar;
    caxis([-7e-4 7e-4])
    ax.FontSize = 30;
    set(gcf, 'Position', [100 100 975 975/1.3333])
    print([print_dir 'amp_sens_' num2str(P(i)) 's'], '-dpng', '-r300')
    close

    phaseK_sens = (A_obs(i)./phase(i)).*H_AK(i,:) + (-B_obs(i)./phase(i)).*H_BK(i,:);
    phaseK_sensmap = reshape(phaseK_sens,num_y,num_x);
    figure
    clf
    ax = gca;
    pc2 = contourf(cgrid{1},cgrid{2}, phaseK_sensmap, 10);
%     set(pc2,'LineStyle','none')
    hold on
    plot(well_locs(:,1), well_locs(:,2), 'ko', 'MarkerFaceColor', 'r', 'MarkerSize', 14)
    xlabel('X Distance (m)')
    ylabel('Y Distance (m)')
    c = colorbar;
    axis equal
    axis([-35 35 -25 45])
    ax.FontSize = 30;
    set(gcf, 'Position', [100 100 975 975/1.3333])
    print([print_dir 'phase_sens_' num2str(P(i)) 's'], '-dpng', '-r300')
    close
end

%% Prior setup - Linear variogram model
distmat_row = dimdist(coords(1,:),coords);
max_dist = max(max((distmat_row(:,:,1).^2 + distmat_row(:,:,2).^2).^.5));
slope = 1/max_dist;

Q_row_corr = slope .*...
            (max_dist - ((distmat_row(:,:,1).^2 + ...
                          distmat_row(:,:,2).^2).^.5));
                      
%% Inversion
R = obs_phasor_covar; R_inv = inv(R);
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
model_err = norm(params_best(1:num_cells) - X(1:num_cells,1)*beta_best(1));
data_err = (y - h_function(params_best))' * R_inv * (y - h_function(params_best));

%% Save output
save_dir = '/.../.../.../';
save_name = ['chkrbrd_' case_name '_10m_wells.mat'];
save([save_dir save_name])