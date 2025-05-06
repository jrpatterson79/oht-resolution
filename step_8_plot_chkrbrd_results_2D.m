% Plot Tomography Results

%% Clean Environment
close all; clear; clc

%% Specify Directory
load_dir = '/Users/jpatt/Dropbox (Personal)/projects/oht_resolution/mat_files/';
print_dir = '/Users/jpatt/Dropbox (Personal)/Apps/Overleaf/patterson_cardiff_2024/fig_comp/';
addpath(genpath('/Users/jpatt/Dropbox (Personal)/projects/oht_resolution/oht/'))

%% Checkerboard Size Analysis
file_name = {'chkrbrd_4_freq_20m_checker';
             'chkrbrd_4_freq';
             'chkrbrd_4_freq_8m_checker';
             };
num_case = numel(file_name);

for w = 1:num_case    
    load([load_dir file_name{w}])
    disp(file_name{w})
    disp(['T_var = ' num2str(lnK_var_est) '; S_var = ' num2str(lnSs_var_est)])
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
    clim([log10(exp(lnK_range(1))) log10(exp(lnK_range(2)))])
end
c = colorbar;
c.Layout.Tile = 'east';
c.Ticks = (-4.2:0.1:-3.8);
c.Label.String = 'log_{10}(T_{est} [m^2/s])';
c.FontSize = 24;
set(gcf, 'Position', [100 100 2025 2025/2.667])
% print([print_dir 'fig_9'], '-djpeg', '-r600')

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
    clim([log10(exp(lnSs_range(1))) log10(exp(lnSs_range(2)))])
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
    disp(['T_var = ' num2str(lnK_var_est) '; S_var = ' num2str(lnSs_var_est)])
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
    % clim([log10(exp(lnK_range(1))) log10(exp(lnK_range(2)))])
    clim([-4.2 -3.8])
end
c = colorbar;
c.Layout.Tile = 'east';
c.Label.String = 'log_{10}(T_{est} [m^2/s])';
c.FontSize = 24;
set(gcf, 'Position', [100 100 2025 2025/2.667])
% print([print_dir 'sing_freq_check'], '-dpng', '-r1200')

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
    clim([log10(exp(lnSs_range(1))) log10(exp(lnSs_range(2)))])
    % title(lbls{j}, 'FontWeight', 'bold', 'FontSize', 30)
end
c = colorbar;
c.Layout.Tile = 'east';
c.Label.String = 'log_{10}(S_{est} [-])';
c.FontSize = 24;
set(gcf, 'Position', [100 100 2025 2025/2.667])
% print([print_dir 'sf_check_s'], '-dpng', '-r300')

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
    disp(['T_var = ' num2str(lnK_var_est) '; S_var = ' num2str(lnSs_var_est)])
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
    % clim([log10(exp(lnK_range(1))) log10(exp(lnK_range(2)))])
end
c = colorbar;
c.Layout.Tile = 'east';
c.Ticks = (-4.2:0.1:-3.8);
c.TickLabels = (-4.2:0.1:-3.8);
c.Label.String = 'log_{10}(T_{est} [m^2/s])';
c.FontSize = 24;
set(gcf, 'Position', [100 100 2025 2025/2.667])
% print([print_dir 'multi_freq_check'], '-djpeg', '-r1200')

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
    clim([log10(exp(lnSs_range(1))) log10(exp(lnSs_range(2)))])
    % title(labels{j}, 'FontWeight', 'bold', 'FontSize', 30)
end
c = colorbar;
c.Layout.Tile = 'east';
c.Label.String = 'log_{10}(S_{est} [-])';
c.FontSize = 24;
set(gcf, 'Position', [100 100 2025 2025/2.667])
% print([print_dir 'mf_check_S'], '-dpng', '-r300')


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
    disp(['T_var = ' num2str(lnK_var_est) '; S_var = ' num2str(lnSs_var_est)])

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
% print([print_dir 'fig_11'], '-djpeg', '-r600')

% Multi-frequency vs 10-m well spacing (Storativity)
figure
clf
tiledlayout(1,num_case, 'TileSpacing', 'loose', 'Padding', 'loose')
for w = 1 : num_case
    nexttile(w)
    load([load_dir file_name{w}])
    disp(file_name{w})
    disp(['T_var = ' num2str(lnK_var_est) '; S_var = ' num2str(lnSs_var_est)])

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
    clim([log10(exp(lnSs_range(1))) log10(exp(lnSs_range(2)))])
end
c = colorbar;
c.Layout.Tile = 'east';
c.Label.String = 'log_{10}(S_{est} [-])';
c.FontSize = 24;
set(gcf, 'Position', [100 100 2025 2025/2.667])
print([print_dir 'fig_S8'], '-dpng', '-r300')
