% Oscillatory Hydraulic Tomography (OHT) Linear Resolution Analysis

% This code loads .mat files created in step_6_inversion_chkrbrd.m and step_7_inversion_chkrbrd_10m_wells.m and then plots the recovered checkerboards for single- and multi-frequency OHT. Error free code execution requires that all single- and multi-freqeuency checkerboard analyses have been conducted and .mat files created. See lines 19, 100, 177, and 256.

% Unmodified, this code will reproduce the Figures 9-11 and Figures S6-S8 in:
% Patterson, J. R., & Cardiff, M. (2025). Multi‐frequency oscillatory hydraulic tomography improves heterogeneity imaging and resolution and reduces uncertainty. Water Resources Research, 61, e2024WR039606. https://doi.org/10.1029/2024WR039606

% This code requires OHT3DINV v.0.16.0 which can be downloaded at https://github.com/wischydro-cardiff/oscillatory-tomography.git

% Code developed by Jeremy Patterson
% Created: March 2022; Updated May 2025

%% Clean Environment
close all; clear; clc

%% Specify Directories
load_dir = '/.../.../.../'; % Directory with inversion output .mat files 
addpath(genpath('/.../.../.../')) % Directory with OHT3DINV

%% Checkerboard Size Analysis
file_name = {'chkrbrd_4_freq_20m_checker';
             'chkrbrd_4_freq';
             'chkrbrd_4_freq_8m_checker';
             };
num_case = numel(file_name);

for w = 1:num_case    
    load([load_dir file_name{w}])
    disp(file_name{w})
    disp(['T_var = ' num2str(lnT_var_est) '; S_var = ' num2str(lnS_var_est)])
    params_plot(:,w) = log10(exp(params_best));
end
% Figure 9 - Variable Checkerboard Size (Transmissivity)
labels = {'20 m checker'; '10 m checker'; '8 m checker'};
subfig_txt = {'A'; 'B'; 'C'};

figure(9)
clf
tiledlayout(1,3, 'TileSpacing', 'loose', 'Padding', 'Loose')
for j = 1 : num_case
    nexttile(j)
    ax = gca;
    p = pcolor(cgrid{1}, cgrid{2}, reshape(params_plot(1:num_cells,j), num_y, num_x));
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
    text(40, 45, subfig_txt{j}, 'FontWeight', 'bold', 'FontSize', 30)
    clim([log10(exp(lnT_range(1))) log10(exp(lnT_range(2)))])
end
c = colorbar;
c.Layout.Tile = 'east';
c.Ticks = (-4.2:0.1:-3.8);
c.Label.String = 'log_{10}(T_{est} [m^2/s])';
c.FontSize = 24;
set(gcf, 'Position', [100 100 2025 2025/2.667])

% Variable Checkerboard Size (Storativity)
subfig_txt = {'A'; 'B'; 'C'};
figure
clf
tiledlayout(1,3, 'TileSpacing', 'loose', 'Padding', 'Loose')
for j = 1 : num_case
    nexttile(j)
    ax = gca;
    p = pcolor(cgrid{1}, cgrid{2}, reshape(params_plot(num_cells+1:2*num_cells,j), num_y, num_x));
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
    text(40, 45, subfig_txt{j}, 'FontWeight', 'bold', 'FontSize', 30)
    clim([log10(exp(lnS_range(1))) log10(exp(lnS_range(2)))])
end
c = colorbar;
c.Layout.Tile = 'east';
c.Label.String = 'log_{10}(S_{est} [-])';
c.FontSize = 24;
set(gcf, 'Position', [100 100 2500 2500/2.667])
% print([print_dir 'fig_S6'], '-dpng', '-r300')

%% Single frequency Checkerboard Analysis
% Load and open files
file_name = {'chkrbrd_1_freq_high';
             'chkrbrd_1_freq_intermed';
             'chkrbrd_1_freq_low';
             };
num_case = numel(file_name);  

for w = 1 : num_case
    load([load_dir file_name{w}])
    disp(file_name{w})
    disp(['T_var = ' num2str(lnT_var_est) '; S_var = ' num2str(lnS_var_est)])
    params_plot(:,w) = log10(exp(params_best));
end

% Figure 9 - Single Frequency Checkerboard (Transmissivity)
lbls = {'P = 10 s'; 'P = 100 s'; 'P = 1600 s'};
figure
clf
tiledlayout(1, num_case, 'TileSpacing', 'loose', 'Padding', 'Loose')
for j = 1 : num_case
    nexttile(j)
    ax = gca;
    p = pcolor(cgrid{1}, cgrid{2}, reshape(params_plot(1:num_cells,j), num_y, num_x));
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
    % clim([log10(exp(lnT_range(1))) log10(exp(lnT_range(2)))])
    clim([-4.2 -3.8])
end
c = colorbar;
c.Layout.Tile = 'east';
c.Label.String = 'log_{10}(T_{est} [m^2/s])';
c.FontSize = 24;
set(gcf, 'Position', [100 100 2025 2025/2.667])

% Figure S5 - Single Frequency Checkerboard (Storativity)
lbls = {'P = 10 s'; 'P = 100 s'; 'P = 1600 s'};
figure
clf
tiledlayout(1, num_case, 'TileSpacing', 'loose', 'Padding', 'Loose')
for j = 1 : num_case
    nexttile(j)
    ax = gca;
    p = pcolor(cgrid{1}, cgrid{2}, reshape(params_plot(num_cells+1:2*num_cells,j), num_y, num_x));
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
    clim([log10(exp(lnS_range(1))) log10(exp(lnS_range(2)))])
    % title(lbls{j}, 'FontWeight', 'bold', 'FontSize', 30)
end
c = colorbar;
c.Layout.Tile = 'east';
c.Label.String = 'log_{10}(S_{est} [-])';
c.FontSize = 24;
set(gcf, 'Position', [100 100 2025 2025/2.667])

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
    disp(['T_var = ' num2str(lnT_var_est) '; S_var = ' num2str(lnS_var_est)])
    params_plot(:,w) = log10(exp(params_best));
end

% Figure 10 - Multi-Frequency checkerboard recovery (Transmissivity)
labels = {'2 Frequencies'; '4 Frequencies'; '6 Frequencies'};
figure(10)
clf
tiledlayout(1, num_case, 'TileSpacing', 'loose', 'Padding', 'Loose')
for j = 1 : num_case
    nexttile(j)
    ax = gca;
    p = pcolor(cgrid{1}, cgrid{2}, reshape(params_plot(1:num_cells,j), num_y, num_x));
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
    clim([-4.2 -3.8])
    % clim([log10(exp(lnT_range(1))) log10(exp(lnT_range(2)))])
end
c = colorbar;
c.Layout.Tile = 'east';
c.Ticks = (-4.2:0.1:-3.8);
c.TickLabels = (-4.2:0.1:-3.8);
c.Label.String = 'log_{10}(T_{est} [m^2/s])';
c.FontSize = 24;
set(gcf, 'Position', [100 100 2025 2025/2.667])

% Figure S6 - Multi-frequency Checkerboard recovery (Storativity)
figure
clf
tiledlayout(1, num_case, 'TileSpacing', 'loose', 'Padding', 'Loose')
for j = 1 : num_case
    nexttile(j)
    ax = gca;
    p = pcolor(cgrid{1}, cgrid{2}, reshape(params_plot(num_cells+1:2*num_cells,j), num_y, num_x));
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
    clim([log10(exp(lnS_range(1))) log10(exp(lnS_range(2)))])
    % title(labels{j}, 'FontWeight', 'bold', 'FontSize', 30)
end
c = colorbar;
c.Layout.Tile = 'east';
c.Label.String = 'log_{10}(S_{est} [-])';
c.FontSize = 24;
set(gcf, 'Position', [100 100 2025 2025/2.667])

%% 10 m well checkerboard analysis
% Load and open files
file_name = {'chkrbrd_6_freq';
             'chkrbrd_1_freq_10m_wells';
             };
num_case = numel(file_name);

% Figure 11 - 10 m well checkerboard analysis (Transmissivity)
titles = {'6 Frequencies'; '1 Frequency'};
subfig_txt = {'A'; 'B'};

figure (11)
clf
tiledlayout(1,num_case, 'TileSpacing', 'loose', 'Padding', 'loose')
for w = 1 : num_case
    nexttile(w)
    load([load_dir file_name{w}])
    disp(file_name{w})
    disp(['T_var = ' num2str(lnT_var_est) '; S_var = ' num2str(lnS_var_est)])

    ax = gca;
    p = pcolor(cgrid{1}, cgrid{2}, reshape(log10(exp(params_best(1:num_cells))), num_y, num_x));
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
    text(45, 47, subfig_txt{w}, 'FontWeight', 'bold', 'FontSize', 30)
    clim([-4.2 -3.8])
end
c = colorbar;
c.Layout.Tile = 'east';
c.Ticks = (-4.2:0.1:-3.8);
c.Label.String = 'log_{10}(T_{est} [m^2/s])';
c.FontSize = 24;
set(gcf, 'Position', [100 100 2025 2025/2.667])

% Multi-frequency vs 10-m well spacing (Storativity)
figure
clf
tiledlayout(1,num_case, 'TileSpacing', 'loose', 'Padding', 'loose')
for w = 1 : num_case
    nexttile(w)
    load([load_dir file_name{w}])
    disp(file_name{w})
    disp(['T_var = ' num2str(lnT_var_est) '; S_var = ' num2str(lnS_var_est)])

    ax = gca;
    p = pcolor(cgrid{1}, cgrid{2}, reshape(log10(exp(params_best(num_cells+1:2*num_cells))), num_y, num_x));
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
    text(45, 47, subfig_txt{w}, 'FontWeight', 'bold', 'FontSize', 30)
    clim([log10(exp(lnS_range(1))) log10(exp(lnS_range(2)))])
end
c = colorbar;
c.Layout.Tile = 'east';
c.Label.String = 'log_{10}(S_{est} [-])';
c.FontSize = 24;
set(gcf, 'Position', [100 100 2025 2025/2.667])