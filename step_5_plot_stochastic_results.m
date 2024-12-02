% Plot Tomography Results

%% Clean Environment
close all; clear; clc

%% Specify Directory
load_dir = '/.../.../'; % Specify directory where output .mat file is located
addpath(genpath('/.../.../')) % Specify directory where OHT function files are located

%% Load File
file_name = {'stoch_1_freq';...
             'stoch_2_freq';...
             'stoch_4_freq';...
             };

num_case = numel(file_name);
for idx = 1 : num_case
    load([load_dir file_name{idx} '.mat'])
    coords_range = (coords(:,1) > -30 & coords(:,1) < 30 & coords(:,2) > -30 & coords(:,2) < 30);
    % Calculate posterior covariance for each realization
    for relz = 1 : num_relz
        param_rmse(relz,idx) = sqrt(mean((log10(exp(params_true_mat(1:num_cells,relz))) - log10(exp(params_best_mat(1:num_cells,relz)))).^2));
        param_square_err(:,relz) = (log10(exp(params_true_mat(1:num_cells,relz))) - log10(exp(params_best_mat(1:num_cells,relz)))).^2;
        param_corr(relz,idx) = corr(params_true_mat(coords_range,relz), params_best_mat(coords_range,relz));
    end
    ensemble_rmse(:,idx) = sqrt(mean(param_square_err,2));
end

%% Stochastic Analysis
% Figure 6 - Parameter correlation / error images per realization
figure (6)
clf
subplot(1,2,1)
ax = gca;
imagesc(param_corr)
ax.XTick = [1;2;3];
ax.XTickLabel = [1;2;4];
xlabel('# of Frequencies')
ylabel('Realization Number')
ax.FontSize = 30;
c = colorbar;
caxis([0.65 1])
c.Label.String = 'Hydraulic Conductivity Correlation';
c.FontSize = 24;

subplot(1,2,2)
ax = gca;
imagesc(param_rmse)
ax.XTick = [1;2;3];
ax.XTickLabel = [1;2;4];
xlabel('# of Frequencies')
% ax.YTickLabel = [];
ax.FontSize = 30;
c = colorbar;
caxis([0 0.8])
c.Label.String = 'Hydraulic Conductivity RMSE (m/s)';
c.FontSize = 24;
set(gcf, 'Position', [100 100 2025/1.1 2025])

% Figure 7 - Ensemble parameter RMSE
figure (7)
clf
tiledlayout(1,3, 'TileSpacing', 'loose', 'Padding', 'Loose')
titles = {'1 Frequency';
          '2 Frequencies';
          '4 Frequencies'};
for i = 1 : num_case
    nexttile(i)
    ax = gca;
    p = pcolor(cgrid{1}, cgrid{2}, reshape(ensemble_rmse(:,i), num_y, num_x));
    p.LineStyle = 'none';
    hold on
    plot(well_locs(:,1), well_locs(:,2), 'ko', 'MarkerSize', 12, 'LineWidth', 2,...
        'MarkerFaceColor', [0.7529 0 0])
    axis equal
    axis([-99 99 -99 99])
    ax.XTick = [-80:40:80]; ax.YTick = ax.XTick;
%     caxis([0 0.8])
    if i == 1
        xlabel('X (m)')
        ylabel('Y (m)')
    end
    ax.FontSize = 30;
    title(titles{i}, 'FontWeight', 'bold', 'FontSize', 30)
end
c = colorbar;
c.Layout.Tile = 'east';
c.Ticks = [0:0.2:0.8];
c.Label.String = 'K_{est} Ensemble Mean RMSE (m/s)';
c.FontSize = 24;
set(gcf,'Position', [100 100 2025 2025/3.1])