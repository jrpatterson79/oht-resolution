% Oscillatory Hydraulic Tomography (OHT) Linear Resolution Analysis

% This code loads .mat files created in step_4_stochastic_inversion.m and plots the inverted tomograms for single- and multi-frequency OHT. Error free code execution requires that all single- and multi-freqeuency stochastic analyses have been conducted and .mat files created. See lines 19 and 110.

% Unmodified, this code will reproduce the Figures 7-8 and Figure S5 in:
% Patterson, J. R., & Cardiff, M. (2025). Multi‐frequency oscillatory hydraulic tomography improves heterogeneity imaging and resolution and reduces uncertainty. Water Resources Research, 61, e2024WR039606. https://doi.org/10.1029/2024WR039606

% Code developed by Jeremy Patterson
% Created: March 2022; Updated May 2025

%% Clean Environment
close all; clear; clc

%% Specify Directories
load_dir = '/.../.../.../'; % Directory with inversion output .mat files 
addpath(genpath('/.../.../.../')) % Directory with OHT3DINV

%% Single Frequency Stochastic Analysis
file_name = {'stoch_1_freq_high';...
             'stoch_1_freq_intermed';...
             'stoch_1_freq_low';...
             };

num_case = numel(file_name);
for idx = 1 : num_case
    load([load_dir file_name{idx} '.mat'])
    coords_range = find(coords(:,1) > -40 & coords(:,1) < 40 & coords(:,2) > -40 & coords(:,2) < 40);
    % coords_range = find(coords(:,1) > -20 & coords(:,1) < 20 & coords(:,2) > -20 & coords(:,2) < 20);

    for relz = 1 : num_relz
        % Transmissivity error metrics
        T_rmse(relz,idx) = sqrt(mean((log10(exp(params_true_mat(1:num_cells,relz))) - log10(exp(params_best_mat(1:num_cells,relz)))).^2));
        T_square_err(:,relz) = (log10(exp(params_true_mat(1:num_cells,relz))) - log10(exp(params_best_mat(1:num_cells,relz)))).^2;
        T_corr(relz,idx) = corr(params_true_mat(coords_range,relz), params_best_mat(coords_range,relz));

        % Storativity error metrics
        S = log10(exp(params_best_mat(num_cells+1:2*num_cells,relz))); S_true = log10(exp(params_true_mat(num_cells+1:2*num_cells,relz)));
        S_rmse(relz,idx) = sqrt(mean((S - S_true).^2));
        S_square_err(:,relz) = (S - S_true).^2;
        S_corr(relz,idx) = corr(params_true_mat(coords_range,relz), S(coords_range));
    end
    % Realizations 16, 17, and 50 are poor quality. Excluded from analysis
    T_square_err = T_square_err(:,[1:15,18:49, 51:55, 57:71, 73:79, 81:92, 94:end]);
    S_square_err = S_square_err(:,[1:15,18:49, 51:55, 57:71, 73:79, 81:92, 94:end]);

    ensemble_T_rmse(:,idx) = sqrt(mean(T_square_err,2));
    ensemble_S_rmse(:,idx) = sqrt(mean(S_square_err,2));
end

% Figure 8 - Single-frequency Transmissivity Ensemble RMSE Maps
figure
clf
tiledlayout(1,3, 'TileSpacing', 'loose', 'Padding', 'Loose')
for i = 1 : num_case
    nexttile(i)
    ax = gca;
    p = pcolor(cgrid{1}, cgrid{2}, reshape(ensemble_T_rmse(:,i), num_y, num_x));
    p.LineStyle = 'none';
    hold on
    plot(well_locs(:,1), well_locs(:,2), 'ko', 'MarkerSize', 12, 'LineWidth', 2,...
        'MarkerFaceColor', [0.7529 0 0])
    axis equal
    axis([-99 99 -99 99])
    ax.XTick = [-80:40:80]; ax.YTick = ax.XTick;
    clim([0 0.8])
    if i == 1
        xlabel('X (m)')
        ylabel('Y (m)')
    end
    ax.FontSize = 30;
end
c = colorbar;
c.Layout.Tile = 'east';
c.Ticks = [0:0.2:0.8];
c.Label.String = 'T_{est} Ensemble Mean RMSE (m^2/s)';
c.FontSize = 24;
set(gcf,'Position', [100 100 2025 2025/3.1])

% Figure S5 - Storativity Ensemble RMSE
figure 
clf
tiledlayout(1,3, 'TileSpacing', 'loose', 'Padding', 'Loose')
for i = 1 : num_case
    nexttile(i)
    ax = gca;
    p = pcolor(cgrid{1}, cgrid{2}, reshape(ensemble_S_rmse(:,i), num_y, num_x));
    p.LineStyle = 'none';
    hold on
    plot(well_locs(:,1), well_locs(:,2), 'ko', 'MarkerSize', 12, 'LineWidth', 2,...
        'MarkerFaceColor', [0.7529 0 0])
    axis equal
    axis([-80 80 -80 80])
    ax.XTick = [-80:40:80]; ax.YTick = ax.XTick;
    clim([0 0.5])
    if i == 1
        xlabel('X (m)')
        ylabel('Y (m)')
    end
    ax.FontSize = 30;
end
c = colorbar;
c.Layout.Tile = 'east';
c.Ticks = [0:0.1:0.5];
c.Label.String = 'S_{est} Ensemble Mean RMSE (-)';
c.FontSize = 24;
set(gcf,'Position', [100 100 2025 2025/3.1])
close all

%% Multi-frequency Stochastic Analysis
file_name = {'stoch_1_freq_low';...
             'stoch_2_freq';...
             'stoch_4_freq';
             'stoch_6_freq'};

num_case = numel(file_name);
for idx = 1 : num_case
    load([load_dir file_name{idx} '.mat'])
    coords_range = find(coords(:,1) > -40 & coords(:,1) < 40 & coords(:,2) > -40 & coords(:,2) < 40);
    % coords_range = find(coords(:,1) > -20 & coords(:,1) < 20 & coords(:,2) > -20 & coords(:,2) < 20);

    for relz = 1 : num_relz
        % Transmissivity error metrics
        T_rmse(relz,idx) = sqrt(mean((log10(exp(params_true_mat(1:num_cells,relz))) - log10(exp(params_best_mat(1:num_cells,relz)))).^2));
        T_square_err(:,relz) = (log10(exp(params_true_mat(1:num_cells,relz))) - log10(exp(params_best_mat(1:num_cells,relz)))).^2;
        T_corr(relz,idx) = corr(params_true_mat(coords_range,relz), params_best_mat(coords_range,relz));

        % Storativity error metrics
        S = log10(exp(params_best_mat(num_cells+1:2*num_cells,relz))); S_true = log10(exp(params_true_mat(num_cells+1:2*num_cells,relz)));
        S_rmse(relz,idx) = sqrt(mean((S - S_true).^2));
        S_square_err(:,relz) = (S - S_true).^2;
        S_corr(relz,idx) = corr(params_true_mat(coords_range,relz), S(coords_range));
    end
    % Realizations 2, 10, 11, 17, and 37 were poor quality inversions -
    % excluded from analysis
    T_square_err = T_square_err(:,[1,3:9,12:16,18:36,38:end]);
    S_square_err = S_square_err(:,[1,3:9,12:16,18:36,38:end]);
    ensemble_T_rmse(:,idx) = sqrt(mean(T_square_err,2));
    ensemble_S_rmse(:,idx) = sqrt(mean(S_square_err,2));
end

% Figure 7 - Transmissivity correlation / error metrics
figure(7)
clf
subplot(1,2,1)
ax = gca;
imagesc(T_corr([1,3:9,12:16,18:36,38:end],:))
ax.XTick = [1;2;3;4];
ax.XTickLabel = [1;2;4;6];
xlabel('# of Frequencies')
ylabel('Realization Number')
ax.FontSize = 30;
c = colorbar;
clim([0.65 1])
c.Label.String = 'Transmissivity Correlation';
c.FontSize = 24;

subplot(1,2,2)
ax = gca;
imagesc(T_rmse([1,3:9,12:16,18:36,38:end],:))
ax.XTick = [1;2;3;4];
ax.XTickLabel = [1;2;4;6];
xlabel('# of Frequencies')
% ax.YTickLabel = [];
ax.FontSize = 30;
c = colorbar;
clim([0 0.8])
c.Label.String = 'Transmissivity RMSE (m^2/s)';
c.FontSize = 24;
set(gcf, 'Position', [100 100 2025/1.1 2025])

%Figure S5 - Storativity error / correlation metrics 
figure
clf
subplot(1,2,1)
ax = gca;
imagesc(S_corr)
ax.XTick = [1;2;3];
ax.XTickLabel = [1;2;4];
xlabel('# of Frequencies')
ylabel('Realization Number')
ax.FontSize = 30;
c = colorbar;
caxis([0.4 1])
c.Label.String = 'Storativity Correlation';
c.FontSize = 24;

subplot(1,2,2)
ax = gca;
imagesc(S_rmse)
ax.XTick = [1;2;3];
ax.XTickLabel = [1;2;4];
xlabel('# of Frequencies')
% ax.YTickLabel = [];
ax.FontSize = 30;
c = colorbar;
caxis([0 0.8])
c.Label.String = 'Storativity RMSE (-)';
c.FontSize = 24;
set(gcf, 'Position', [100 100 2025/1.1 2025])

% Figure 8 - Multi-frequency Transmissivity Ensemble RMSE
figure
clf
tiledlayout(1,3, 'TileSpacing', 'loose', 'Padding', 'Loose')
for i = 2 : num_case
    nexttile(i-1)
    ax = gca;
    p = pcolor(cgrid{1}, cgrid{2}, reshape(ensemble_T_rmse(:,i), num_y, num_x));
    p.LineStyle = 'none';
    hold on
    plot(well_locs(:,1), well_locs(:,2), 'ko', 'MarkerSize', 12, 'LineWidth', 2,...
        'MarkerFaceColor', [0.7529 0 0])
    axis equal
    axis([-99 99 -99 99])
    ax.XTick = [-80:40:80]; ax.YTick = ax.XTick;
    caxis([0 0.8])
    if i == 2
        xlabel('X (m)')
        ylabel('Y (m)')
    end
    ax.FontSize = 30;
end
c = colorbar;
c.Layout.Tile = 'east';
c.Ticks = [0:0.2:0.8];
c.Label.String = 'T_{est} Ensemble Mean RMSE (m^2/s)';
c.FontSize = 24;
set(gcf,'Position', [100 100 2025 2025/3.1])

% Figure S5 - Multi-Frequency Storativity Ensemble RMSE
figure 
clf
tiledlayout(1,3, 'TileSpacing', 'loose', 'Padding', 'Loose')
for i = 2 : num_case
    nexttile(i-1)
    ax = gca;
    p = pcolor(cgrid{1}, cgrid{2}, reshape(ensemble_S_rmse(:,i), num_y, num_x));
    p.LineStyle = 'none';
    hold on
    plot(well_locs(:,1), well_locs(:,2), 'ko', 'MarkerSize', 12, 'LineWidth', 2,...
        'MarkerFaceColor', [0.7529 0 0])
    axis equal
    axis([-80 80 -80 80])
    ax.XTick = [-80:40:80]; ax.YTick = ax.XTick;
    clim([0 0.5])
    if i == 2
        xlabel('X (m)')
        ylabel('Y (m)')
    end
    ax.FontSize = 30;
end
c = colorbar;
c.Layout.Tile = 'east';
c.Ticks = [0:0.1:0.5];
c.Label.String = 'S_{est} Ensemble Mean RMSE (-)';
c.FontSize = 24;
set(gcf,'Position', [100 100 2025 2025/3.1])