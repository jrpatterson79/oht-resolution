% Plot Checkerboard Results

%% Clean Environment
close all; clear; clc

%% Specify Directory
load_dir = '/.../.../'; % Directory where output .mat files are located
addpath(genpath('/.../.../')) % Directory where OHT function files are located

%% Checkerboard Size Analysis
file_name = {'chkrbrd_4_freq_20m_checker';
             'chkrbrd_4_freq';
             'chkrbrd_4_freq_8m_checker';
             };
num_case = numel(file_name);

for w = 1:num_case
    load([load_dir file_name{w}])
    disp(file_name{w})
    disp(['K_var = ' num2str(lnK_var_est) '; Ss_var = ' num2str(lnSs_var_est)])
    % Calculate Pearson correlation coefficient
    param_corr(w,:) = corr(params_true(coords_range), params_best(coords_range));
    params_plot(:,w) = params_best;
end

% Figure 8 - Variable Checkerboard Size
labels = {'20 m checker'; '10 m checker'; '8 m checker'};
figure(8)
clf
tiledlayout(1,3, 'TileSpacing', 'loose', 'Padding', 'Loose')
for j = 1 : num_case
    nexttile(j)
    ax = gca;
    p = pcolor(cgrid{1}, cgrid{2}, reshape(log10(exp(params_plot(1:num_cells,j))), num_y, num_x));
    hold on
    plot(well_locs(:,1), well_locs(:,2), 'ko', 'LineWidth', 2,...
        'MarkerFaceColor', [0.7529 0 0], 'MarkerSize', 12)
    p.LineStyle = 'none';
    axis equal
    axis([-50 50 -50 50])
    ax.XTick = [-50:25:50]; ax.YTick = [-50:25:50];
    if j == 1
        xlabel('X (m)')
        ylabel('Y (m)')
    end
    ax.FontSize = 30;
    title(labels{j}, 'FontWeight', 'bold', 'FontSize', 30)
    caxis([log10(exp(lnK_range(1))) log10(exp(lnK_range(2)))])
end
c = colorbar;
c.Layout.Tile = 'east';
c.Label.String = 'log_{10}(K_{est} [m/s])';
c.FontSize = 24;
set(gcf, 'Position', [100 100 2500 2500/2.667])

%% Single frequency Checkerboard Analysis
% Load and open files
file_name = {'chkrbrd_1_freq_high';
             'chkrbrd_1_freq';
             'chkrbrd_1_freq_low';
             };
num_case = numel(file_name);  

for w = 1 : num_case
    load([load_dir file_name{w}])
    disp(file_name{w})
    disp(['K_var = ' num2str(lnK_var_est) '; Ss_var = ' num2str(lnSs_var_est)])
    % Calculate Pearson correlation coefficient 
    coords_range = (coords(:,1) > -30 & coords(:,1) < 30 & coords(:,2) > -30 & coords(:,2) < 30);
    param_rmse_sing_freq(w,:) = sqrt(mean((params_best(coords_range)-params_true(coords_range)).^2));
    params_plot(:,w) = params_best;
end

% Figure 9 - Single Frequency Checkerboard
lbls = {'P = 10 s'; 'P = 100 s'; 'P = 1600 s'};
figure(9)
clf
tiledlayout(1, num_case, 'TileSpacing', 'loose', 'Padding', 'Loose')
for j = 1 : num_case
    nexttile(j)
    ax = gca;
    p = pcolor(cgrid{1}, cgrid{2}, reshape(log10(exp(params_plot(1:num_cells,j))), num_y, num_x));
    hold on
    plot(well_locs(:,1), well_locs(:,2), 'ko', 'LineWidth', 2,...
        'MarkerFaceColor', [0.7529 0 0], 'MarkerSize', 12)
    p.LineStyle = 'none';
    axis equal
    axis([-50 50 -50 50])
    ax.XTick = [-50:25:50]; ax.YTick = [-50:25:50];
    if j == 1
        xlabel('X (m)')
        ylabel('Y (m)')
    end
    ax.FontSize = 30;
    caxis([log10(exp(lnK_range(1))) log10(exp(lnK_range(2)))])
    title(lbls{j}, 'FontWeight', 'bold', 'FontSize', 30)
end
c = colorbar;
c.Layout.Tile = 'east';
c.Label.String = 'log_{10}(K_{est} [m/s])';
c.FontSize = 24;
set(gcf, 'Position', [100 100 2500 2500/2.667])

%% Multi-frequency Checkerboard Analysis
% Load and open files
file_name = {'chkrbrd_2_freq';
             'chkrbrd_4_freq';
             'chkrbrd_6_freq';
             };
num_case = numel(file_name);  

for w = 1:num_case
    load([load_dir file_name{w}])
    disp(file_name{w})
    disp(['K_var = ' num2str(lnK_var_est) '; Ss_var = ' num2str(lnSs_var_est)])
    % Calculate Pearson correlation coefficient
    coords_range = (coords(:,1) > -30 & coords(:,1) < 30 & coords(:,2) > -30 & coords(:,2) < 30);
    param_rmse_multi_freq(w,:) = sqrt(mean((params_best(coords_range)-params_true(coords_range)).^2));    
    params_plot(:,w) = params_best;
end

% Figure 10 - Multi-frequency Checkerboard recovery
labels = {'2 Frequencies'; '4 Frequencies'; '6 Frequencies'};
figure(10)
clf
tiledlayout(1, num_case, 'TileSpacing', 'loose', 'Padding', 'Loose')
for j = 1 : num_case
    nexttile(j)
    ax = gca;
    p = pcolor(cgrid{1}, cgrid{2}, reshape(log10(exp(params_plot(1:num_cells,j))), num_y, num_x));
    hold on
    plot(well_locs(:,1), well_locs(:,2), 'ko', 'LineWidth', 2,...
        'MarkerFaceColor', [0.7529 0 0], 'MarkerSize', 12)
    p.LineStyle = 'none';
    axis equal
    axis([-50 50 -50 50])
    ax.XTick = [-50:25:50]; ax.YTick = [-50:25:50];
    if j == 1
        xlabel('X (m)')
        ylabel('Y (m)')
    end
    ax.FontSize = 30;
    caxis([log10(exp(lnK_range(1))) log10(exp(lnK_range(2)))])
    title(labels{j}, 'FontWeight', 'bold', 'FontSize', 30)
end
c = colorbar;
c.Layout.Tile = 'east';
c.Label.String = 'log_{10}(K_{est} [m/s])';
c.FontSize = 24;
set(gcf, 'Position', [100 100 2500 2500/2.667])

%% 10 m well checkerboard analysis
% Load and open files
file_name = {'chkrbrd_6_freq';
             'chkrbrd_1_freq_10m_wells';
             };
num_case = numel(file_name);

% Figure 11 - 10 m well checkerboard analysis
titles = {'6 Frequencies'; '1 Frequency'};
figure (11)
clf
tiledlayout(1,num_case, 'TileSpacing', 'loose', 'Padding', 'loose')
for w = 1:num_case
    load([load_dir file_name{w}])
    disp(file_name{w})
    disp(['K_var = ' num2str(lnK_var_est) '; Ss_var = ' num2str(lnSs_var_est)])
    % Calculate Pearson correlation coefficient
    param_corr(w,:) = corr(params_true(coords_range), params_best(coords_range));
    params_plot(:,w) = params_best;
    
    nexttile(w)
    ax = gca;
    p = pcolor(cgrid{1}, cgrid{2}, reshape(log10(exp(params_plot(1:num_cells,w))), num_y, num_x));
    hold on
    plot(well_locs(:,1), well_locs(:,2), 'ko', 'LineWidth', 2,...
        'MarkerFaceColor', [0.7529 0 0], 'MarkerSize', 12)
    p.LineStyle = 'none';
    axis equal
    axis([-50 50 -50 50])
    ax.XTick = [-50:25:50]; ax.YTick = [-50:25:50];
    if w == 1
        xlabel('X (m)')
        ylabel('Y (m)')
    end
    ax.FontSize = 30;
    title(titles{w}, 'FontSize', 30, 'FontWeight', 'bold')
    caxis([log10(exp(lnK_range(1))) log10(exp(lnK_range(2)))])
end
c = colorbar;
c.Layout.Tile = 'east';
c.Label.String = 'log_{10}(K_{est} [m/s])';
c.FontSize = 24;
set(gcf, 'Position', [100 100 2025 2025/2.667])