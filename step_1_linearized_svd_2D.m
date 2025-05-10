% Oscillatory Hydraulic Tomography Linear Resolution Analysis

% This code conducts numerical modeling experiments to assess transmissivity and storativity resolution associated with OHT using linear inverse modeling theory. This code numerically approximates the Jacobian matrix with homogeneous aquifer parameters then conducts Singular Value Decomposition to plot the singular value spectra and resolution matrices for single- and multi-frequency OHT testing.

% This code generates Figure 1 and Figure 2 as seen in:
% Patterson, J. R., & Cardiff, M. (2025). Multi‐frequency oscillatory hydraulic tomography improves heterogeneity imaging and resolution and reduces uncertainty. Water Resources Research, 61, e2024WR039606. https://doi.org/10.1029/2024WR039606

% This code requires OHT3DINV v.0.16.0 which can be downloaded at https://github.com/wischydro-cardiff/oscillatory-tomography

% Code developed by Jeremy Patterson
% Created: March 2022; Updated May 2025

%% Clean Environment
clear all; close all; clc;

% Add path to OHT v0.0.16
addpath(genpath('/.../.../.../'))

%% Initialize Forward Modeling Setup
%PARAMETERS: Visualization and saving options
pause_length = 0.5;

%% Specify Model Domain Geometry
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
bdry_L = zeros(6,1);
bdrys = struct('types',bdry_types,'vals',bdry_vals,'leaks',bdry_L);

% Locations of all wells (pumping and observation)
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

%% Define OHT Testing Parameters
% List defining each observation. Columns are:
% pumping angular frequency, pumping well, Q_max, and observation well
% NOTE: The model will run the fastest if this list is sorted by angular
% frequency and pumping well (the first two columns)

V_total = 0.01;
P = [10 50 90 300 700 1600];
test_list = [];
for i = 1:numel(P)
    for j = 1:num_wells
        for k = j+1:num_wells
            test_list = [...
                test_list; ...
                (2*pi)/P(i) j V_total*pi/P(i) k ...
                ];
        end
    end
end

%% Create the "official" input files needed by the steady periodic model.

%This step automatically creates the input files needed to run all of
%the steady-periodic models.

[experiment] = OHT_create_inputs(well_locs,test_list,domain);

% Calculate number of observations, size of grid, and set up a
% meshgridded set of points to do plotting of results.
num_omegas = size(experiment,1);
num_obs = size(test_list,1);
num_x = numel(domain.x) - 1;
num_y = numel(domain.y) - 1;
num_cells = num_x*num_y;

%% Setup grid of cell centers for plotting and geostatistical setups
[coords, cgrid] = plaid_cellcenter_coord(domain);
distmat_row = dimdist(coords(1,:),coords); % Distance matrix for correlation matrix

fig_srcobs_weights = figure;
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
        hold on
        plot(well_locs(:,1), well_locs(:,2), 'ko', 'MarkerFaceColor', [0.7592 0 0], 'MarkerSize', 12)
        colorbar
        title(['Period group ', num2str(i), ', Pump weights, test/obs pair #', num2str(j)]);
        axis equal
        axis([-40 40 -40 40])

        subplot(1,2,2)
        obsmap = reshape(full(experiment(i).obs(:,obs_loc)),num_y,num_x);
        om = pcolor(cgrid{1},cgrid{2},obsmap);
        set(om,'LineStyle','none')
        hold on
        plot(well_locs(:,1), well_locs(:,2), 'ko', 'MarkerFaceColor', [0.7592 0 0], 'MarkerSize', 12)
        colorbar
        title(['Period group ', num2str(i), ', Observation weights, test/obs pair #', num2str(j)]);
        axis equal
        axis([-40 40 -40 40])
        pause(pause_length)
    end
end

%% True parameter field for linearized analysis
% PARAMETERS: Associated with creating test case for fields
%Used by both methods
lnT_mean = -9.2;
lnS_mean = -11.2;
params_true = [lnT_mean*ones(num_cells,1); lnS_mean*ones(num_cells,1)];

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
y_func = @(params) OHT_run_distribKSs(params, domain, bdrys, experiment, 1);
H_tilde_func = @(params) OHT_run_distribKSs(params, ...
    domain,bdrys,experiment,3);

%Generate example sensitivity map given initial homogeneous parameter
%guesses
tic
H_tilde = H_tilde_func(params_true);
toc

%%  Singular Value Decomposition
for i = 1:num_omegas
    if i ~= 1
        labels{i} = [num2str(i), ' Frequencies'];
    else
        labels{i} = [num2str(i), ' Frequency'];
    end
    % Parse full Jacobian based on the number of frequency components and
    % conduct SVD
    H_partial = H_tilde(1:72*i,:);
    [U,S,V] = svd(H_partial);

    % Extract singular values (Normalized by largest singular value)
    sing_val_mat(1:72*i,i) = diag(S) ./ max(diag(S)); 

    % Apply singular value threshold for reduced rank approximation
    p = 72*i;
    Vp = V(:,1:p);
    R = Vp * Vp';
    R_mat(:,i) = diag(R);

end

%% Figures
% Figure 1 - Singular value decomposition
figure(1)
clf 
ax = gca;
plot(sing_val_mat, 'LineWidth', 2);
ax.YScale = 'log';
grid on
xlabel('Singular Value Index')
ylabel('Singular Value')
axis([0 500 1e-15 1e0])
ax.FontSize = 30;
ax.MinorGridAlpha = 0.7;
l = legend(labels);
l.Location = 'Southwest';
l.FontSize = 24;
set(gcf, 'Position', [100 100 975 975/1.3333])

% Figure 1 inset
figure
clf 
ax = gca;
plot(sing_val_mat, 'LineWidth', 2);
ax.YScale = 'log';
grid on
xlabel('Singular Value Index')
ylabel('Singular Value')
axis([0 100 1e-5 1e0])
ax.FontSize = 30;
ax.MinorGridAlpha = 0.7;
set(gcf, 'Position', [100 100 975 975/1.3333])

% Figure 2 - Transmissivity Resolution matrix
figure(2)
clf
t = tiledlayout(2,3);
for j = 1 : numel(P)
    ax = nexttile;
    p = pcolor(cgrid{1}, cgrid{2}, reshape(R_mat(1:num_cells,j),num_y,num_x));
    p.LineStyle = 'none';
    hold on
    plot(well_locs(:,1), well_locs(:,2), 'ko', 'MarkerSize', 12, 'MarkerFaceColor', [0.7592 0 0])
    axis([-40 40 -40 40])
    if j == 1 || j == 4
        xlabel('X (m)')
        ylabel('Y (m)')
    end
    ax.FontSize = 30;
    caxis([0 0.4])
    title(labels{j}, 'FontWeight', 'bold', 'FontSize', 30)
end
t.TileSpacing = 'loose';
t.Padding = 'loose';
c = colorbar;
c.Layout.Tile = 'east';
c.Label.String = 'Resolution (-)';
c.FontSize = 24;
set(gcf, 'Position', [100 100 2025 2025/1.75])

% Figure S1 - Storativity Resolution matrix
figure
clf
t = tiledlayout(2,3);
for j = 1 : numel(P)
    ax = nexttile;
    p = pcolor(cgrid{1}, cgrid{2}, reshape(R_mat(num_cells+1:end,j),num_y,num_x));
    p.LineStyle = 'none';
    hold on
    plot(well_locs(:,1), well_locs(:,2), 'ko', 'MarkerSize', 12, 'MarkerFaceColor', [0.7592 0 0])
    axis([-40 40 -40 40])
    if j == 1 || j == 4
        xlabel('X (m)')
        ylabel('Y (m)')
    end
    ax.FontSize = 30;
    caxis([0 0.6])
    title(labels{j}, 'FontWeight', 'bold', 'FontSize', 30)
end
t.TileSpacing = 'loose';
t.Padding = 'loose';
c = colorbar;
c.Layout.Tile = 'east';
c.Label.String = 'Resolution (-)';
c.FontSize = 24;
set(gcf, 'Position', [100 100 2025 2025/1.75])