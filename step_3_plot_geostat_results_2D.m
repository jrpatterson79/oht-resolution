% Plot Tomography Results

%% Clean Environment
close all; clear; clc

%% Specify Directory
load_dir = '/.../.../.../'; % Directory with inversion output .mat files 
addpath(genpath('/.../.../.../'))

%% Single-Frequency Analysis
file_name = {'geostat_1_freq_high';
             'geostat_1_freq_intermed';
             'geostat_1_freq_low'};

num_case = numel(file_name);

for w = 1 : num_case
    % Load the inversion results
    reg_idx = 9;
    load([load_dir file_name{w}])
    disp(file_name{w})
    disp(['T_var = ' num2str(lnK_var_est(reg_idx)) '; S_var = ' num2str(lnSs_var_est(reg_idx))])
    
    % Parse T and S parameters for plotting - Calculate correlation
    % coefficient
    T_range = log10(exp((lnK_range))); S_range = log10(exp(lnSs_range));
    coords_range = find(coords(:,1) > -40 & coords(:,1) < 40 & coords(:,2) > -40 & coords(:,2) < 40);
    % coords_range = find(coords(:,1) > -20 & coords(:,1) < 20 & coords(:,2) > -20 & coords(:,2) < 20);

    T_est_plot(:,w) = log10(exp(params_best_mat(1:num_cells,reg_idx)));
    S_est_plot(:,w) = log10(exp(params_best_mat(num_cells+1:2*num_cells,reg_idx)));
    S_true = log10(exp(params_true(num_cells+1:2*num_cells)));

    param_corr(w,:) = [corr(params_true(coords_range), T_est_plot(coords_range, w))
                       corr(S_true(coords_range), S_est_plot(coords_range, w))];
    
    

    %% Calculate Posterior Covariance (T and S) 
    reg_param = [lnK_var_est(reg_idx); lnSs_var_est(reg_idx)];
    num_drift = size(X,2);

    H_local = H_tilde_func(params_best_mat(:,reg_idx));    
    QHt = zeros(2*num_cells, 2*num_obs);
    Q_row_reg = Q_row_corr .* reg_param;
    Qprod_func = @(invec) covar_product_K_Ss(Q_row_corr.*reg_param(1), Q_row_corr.*reg_param(2),invec,num_x,num_y);

    for k = 1 : 2*num_obs
        QHt(:,k) = Qprod_func(H_local(k,:)');
    end

    % This calculates the posterior covariance for K only
    HQHt = [H_local*QHt + R];
    HX = H_local * X;
    QHtX = [QHt X];

    covmult_part = [HQHt HX; HX' zeros(num_drift,num_drift)] \ QHtX';
    prior_var = [Q_row_reg(1,1) .* ones(num_cells,1); Q_row_reg(2,1) .* ones(num_cells,1)];

    for i = 1 : 2*num_cells
        post_var(i,w) = log10(exp(prior_var(i) - QHtX(i,:)*covmult_part(:,i)));
    end 
end

%% Transmissivity Figures
% Figure 4 - Single-frequency transmissivity tomogram with parameter cross-plot
for j =  1 : num_case
    figure
    clf
    subplot(1,2,1)
    ax = gca;
    p = pcolor(cgrid{1}, cgrid{2}, reshape(T_est_plot(:,j), num_y, num_x));
    hold on
    plot(well_locs(:,1), well_locs(:,2), 'ko', 'LineWidth', 2,...
        'MarkerFaceColor', [0.7529 0 0], 'MarkerSize', 12)
    p.LineStyle = 'none';
    hold off
    axis equal
    axis([-50 50 -50 50])
    ax.XTick = (-50:25:50); ax.YTick = [-50:25:50];
    xlabel('X (m)')
    ylabel('Y (m)')
    ax.FontSize = 30;
    clim([T_range(1) T_range(2)])
    c = colorbar;
    c.Label.String = 'log_{10}(T_{est} [m^2/s])';
    c.FontSize = 28;

    subplot(1,2,2)
    ax = gca;
    plot([-6 -2], [-6 -2], '-', 'Color', 0.5*[1 1 1], 'LineWidth', 3)
    hold on
    plot(log10(exp(params_true(coords_range))), T_est_plot(coords_range,j), 'ko', ...
        'MarkerFaceColor', [0 0.4470 0.7410], 'MarkerSize', 12)
    axis equal
    axis([-6 -2 -6 -2])
    ax.YTick = [-6:1:-2];
    xlabel('log_{10}(T_{true} [m^2/s])')
    ylabel('log_{10}(T_{est} [m^2/s])')
    ax.FontSize = 30;
    ax.YAxisLocation = 'right';
    text(-5.8, -2.3, ['\rho = ' num2str(round(param_corr(j,1), 2))], 'FontSize', 30, 'FontWeight', 'bold')
    set(gcf, 'Position', [100 100 2025 2025/2.6667])
end

% Figure 6 - Transmissivity posterior covariance (uncertainty)
figure(6)
clf 
t = tiledlayout(1,3);
for j = 1 : num_case
    ax = nexttile;
    p = pcolor(cgrid{1}, cgrid{2}, reshape(sqrt(post_var(1:num_cells,j)), num_y, num_x));
    p.LineStyle = 'none';
    hold on
    plot(well_locs(:,1), well_locs(:,2), 'ko', 'LineWidth', 2,...
        'MarkerSize', 12, 'MarkerFaceColor', [0.7529 0 0])
    hold off
    axis equal
    axis square
    axis([-95 95 -95 95])
    ax.XTick = [-80:40:80]; ax.YTick = ax.XTick;
%     c = colorbar;
    caxis([0.04 0.25])
    if j == 1
        xlabel('X (m)')
        ylabel('Y (m)')
    end
    ax.FontSize = 30;
    % title(labels{j}, 'FontWeight', 'bold', 'FontSize', 30)
end
t.TileSpacing = 'loose';
t.Padding = 'loose';
c = colorbar;
c.Layout.Tile = 'east';
c.Label.String = '\sigma log_{10}(T [m^2/s])';
c.FontSize = 24;
set(gcf, 'Position', [100 100 2025 2025/2.6667])


% Figure S2 - Single-frequency storativity tomograms with parameter cross-plot
for j =  1 : num_case
    figure
    clf
    %Storativity tomogram
    subplot(1,2,1)
    ax = gca;
    p = pcolor(cgrid{1}, cgrid{2}, reshape(S_est_plot(:,j), num_y, num_x));
    hold on
    plot(well_locs(:,1), well_locs(:,2), 'ko', 'LineWidth', 2,...
        'MarkerFaceColor', [0.7529 0 0], 'MarkerSize', 12)
    p.LineStyle = 'none';
    hold off
    axis equal
    axis([-50 50 -50 50])
    ax.XTick = [-50:25:50]; ax.YTick = [-50:25:50];
    xlabel('X (m)')
    ylabel('Y (m)')
    ax.FontSize = 30;
    caxis([S_range(1) S_range(2)])
    c = colorbar;
    c.Label.String = 'log_{10}(S_{est} [-])';
    c.FontSize = 28;
    
    % Figure S2 - Storativity cross-plot
    subplot(1,2,2)
    ax = gca;
    plot([-6 -4], [-6 -4], '-', 'Color', 0.5*[1 1 1], 'LineWidth', 3)
    hold on
    plot(S_true(coords_range), S_est_plot(coords_range,j), 'ko', ...
        'MarkerFaceColor', [0 0.4470 0.7410], 'MarkerSize', 12)
    axis equal
    axis([-6 -4 -6 -4])
    ax.YTick = [-6:0.5:-4];
    xlabel('log_{10}(S_{true} [m/s])')
    ylabel('log_{10}(S_{est} [m/s])')
    ax.FontSize = 30;
    ax.YAxisLocation = 'right';
    text(-5.9, -4.15, ['\rho = ' num2str(round(param_corr(j,2), 2))], 'FontSize', 30, 'FontWeight', 'bold')
    set(gcf, 'Position', [100 100 2025 2025/2.6667])
end

% Figure S2 - Storativity posterior covariance
figure
clf 
t = tiledlayout(1,3);
for j = 1 : num_case
    ax = nexttile;
    p = pcolor(cgrid{1}, cgrid{2}, reshape(sqrt(post_var(num_cells+1:2*num_cells,j)), num_y, num_x));
    p.LineStyle = 'none';
    hold on
    plot(well_locs(:,1), well_locs(:,2), 'ko', 'LineWidth', 2,...
        'MarkerSize', 12, 'MarkerFaceColor', [0.7529 0 0])
    hold off
    axis equal
    axis square
    axis([-95 95 -95 95])
    ax.XTick = [-80:40:80]; ax.YTick = ax.XTick;
    clim([0.05 0.3])
    if j == 1
        xlabel('X (m)')
        ylabel('Y (m)')
    end
    ax.FontSize = 30;
    % title(labels{j}, 'FontWeight', 'bold', 'FontSize', 30)
end
t.TileSpacing = 'loose';
t.Padding = 'loose';
c = colorbar;
c.Layout.Tile = 'east';
c.Label.String = '\sigma log_{10}(S [-])';
c.FontSize = 24;
set(gcf, 'Position', [100 100 2025 2025/2.6667])
% print([print_dir 'fig_S4B'], '-dpng', '-r300')


%% Multi-Frequency Analysis
file_name = {'geostat_2_freq';
             'geostat_4_freq';
             'geostat_6_freq';
             };
num_case = numel(file_name);

for w = 1 : num_case
    % Load the inversion results
    reg_idx = 9;
    load([load_dir file_name{w}])
    disp(file_name{w})
    disp(['T_var = ' num2str(lnK_var_est(reg_idx)) '; S_var = ' num2str(lnSs_var_est(reg_idx))])
    
    % Parse T and S parameters for plotting - Calculate correlation
    % coefficient
    T_range = log10(exp((lnK_range))); S_range = log10(exp(lnSs_range));
    coords_range = find(coords(:,1) > -40 & coords(:,1) < 40 & coords(:,2) > -40 & coords(:,2) < 40);
    % coords_range = find(coords(:,1) > -20 & coords(:,1) < 20 & coords(:,2) > -20 & coords(:,2) < 20);

    T_est_plot(:,w) = log10(exp(params_best_mat(1:num_cells,reg_idx)));
    S_est_plot(:,w) = log10(exp(params_best_mat(num_cells+1:2*num_cells,reg_idx)));
    S_true = log10(exp(params_true(num_cells+1:2*num_cells)));

    param_corr(w,:) = [corr(params_true(coords_range), T_est_plot(coords_range, w))
                       corr(S_true(coords_range), S_est_plot(coords_range, w))];
    
    

    %% Calculate Posterior Covariance (T and S) 
    reg_param = [lnK_var_est(reg_idx); lnSs_var_est(reg_idx)];
    num_drift = size(X,2);

    H_local = H_tilde_func(params_best_mat(:,reg_idx));    
    QHt = zeros(2*num_cells, 2*num_obs);
    Q_row_reg = Q_row_corr .* reg_param;
    Qprod_func = @(invec) covar_product_K_Ss(Q_row_corr.*reg_param(1), Q_row_corr.*reg_param(2),invec,num_x,num_y);

    for k = 1 : 2*num_obs
        QHt(:,k) = Qprod_func(H_local(k,:)');
    end

    HQHt = [H_local*QHt + R];
    HX = H_local * X;
    QHtX = [QHt X];

    covmult_part = [HQHt HX; HX' zeros(num_drift,num_drift)] \ QHtX';
    prior_var = [Q_row_reg(1,1) .* ones(num_cells,1); Q_row_reg(2,1) .* ones(num_cells,1)];

    for i = 1 : 2*num_cells
        post_var(i,w) = log10(exp(prior_var(i) - QHtX(i,:)*covmult_part(:,i)));
    end 
    norm_err(:,w) = [(T_est_plot(:,w) - log10(exp(params_true(1:num_cells)))) ./ sqrt(post_var(1:num_cells,w));
                     (S_est_plot(:,w) - log10(exp(params_true(num_cells+1:2*num_cells)))) ./ sqrt(post_var(num_cells+1:2*num_cells,w))];
end

% Figure 5 - Multi-frequency transmissivity tomogram with parameter cross-plot
labels = {'2 Frequency'; '4 Frequencies'; '6 Frequencies'};
for j =  1 : num_case
    figure
    clf
    subplot(1,2,1)
    ax = gca;
    p = pcolor(cgrid{1}, cgrid{2}, reshape(T_est_plot(:,j), num_y, num_x));
    hold on
    plot(well_locs(:,1), well_locs(:,2), 'ko', 'LineWidth', 2,...
        'MarkerFaceColor', [0.7529 0 0], 'MarkerSize', 12)
    p.LineStyle = 'none';
    hold off
    axis equal
    axis([-50 50 -50 50])
    ax.XTick = [-50:25:50]; ax.YTick = [-50:25:50];
    xlabel('X (m)')
    ylabel('Y (m)')
    ax.FontSize = 30;
    clim([T_range(1) T_range(2)])
    c = colorbar;
    c.Label.String = 'log_{10}(T_{est} [m^2/s])';
    c.FontSize = 28;

    subplot(1,2,2)
    ax = gca;
    plot([-6 -2], [-6 -2], '-', 'Color', 0.5*[1 1 1], 'LineWidth', 3)
    hold on
    plot(log10(exp(params_true(coords_range))), T_est_plot(coords_range,j), 'ko', ...
        'MarkerFaceColor', [0 0.4470 0.7410], 'MarkerSize', 12)
    axis equal
    axis([-6 -2 -6 -2])
    ax.YTick = [-6:1:-2];
    xlabel('log_{10}(T_{true} [m^2/s])')
    ylabel('log_{10}(T_{est} [m^2/s])')
    ax.FontSize = 30;
    ax.YAxisLocation = 'right';
    text(-5.8, -2.3, ['\rho = ' num2str(round(param_corr(j,1), 2))], 'FontSize', 30, 'FontWeight', 'bold')
    set(gcf, 'Position', [100 100 2025 2025/2.6667])
end

% Figure 6 - Multi-frequency transmissivity posterior covariance
figure(6)
clf 
t = tiledlayout(1,3);
for j = 1 : num_case
    ax = nexttile;
    p = pcolor(cgrid{1}, cgrid{2}, reshape(sqrt(post_var(1:num_cells,j)), num_y, num_x));
    p.LineStyle = 'none';
    hold on
    plot(well_locs(:,1), well_locs(:,2), 'ko', 'LineWidth', 2,...
        'MarkerSize', 12, 'MarkerFaceColor', [0.7529 0 0])
    hold off
    axis equal
    axis square
    axis([-95 95 -95 95])
    ax.XTick = [-80:40:80]; ax.YTick = ax.XTick;
    clim([0.04 0.25])
    if j == 1
        xlabel('X (m)')
        ylabel('Y (m)')
    end
    ax.FontSize = 30;
    % title(labels{j}, 'FontWeight', 'bold', 'FontSize', 30)
end
t.TileSpacing = 'loose';
t.Padding = 'loose';
c = colorbar;
% c.Ticks = [0.03 0.09 0.16 0.23 0.3];
c.Layout.Tile = 'east';
c.Label.String = '\sigma log_{10}(T [m^2/s])';
c.FontSize = 24;
set(gcf, 'Position', [100 100 2025 2025/2.6667])

% Figure S3 - Storativity tomogram with parameter cross-plot
labels = {'2 Frequency'; '4 Frequencies'; '6 Frequencies'};
for j =  1 : num_case
    figure
    clf

    % Figure S3 - Multi-frequency storativity tomogram
    subplot(1,2,1)
    ax = gca;
    p = pcolor(cgrid{1}, cgrid{2}, reshape(S_est_plot(:,j), num_y, num_x));
    hold on
    plot(well_locs(:,1), well_locs(:,2), 'ko', 'LineWidth', 2,...
        'MarkerFaceColor', [0.7529 0 0], 'MarkerSize', 12)
    p.LineStyle = 'none';
    hold off
    axis equal
    axis([-50 50 -50 50])
    ax.XTick = [-50:25:50]; ax.YTick = [-50:25:50];
    xlabel('X (m)')
    ylabel('Y (m)')
    ax.FontSize = 30;
    clim([S_range(1) S_range(2)])
    c = colorbar;
    c.Label.String = 'log_{10}(S_{est} [-])';
    c.FontSize = 28;
    
    % Storativity cross-plot
    subplot(1,2,2)
    ax = gca;
    plot([-6 -4], [-6 -4], '-', 'Color', 0.5*[1 1 1], 'LineWidth', 3)
    hold on
    plot(S_true(coords_range), S_est_plot(coords_range,j), 'ko', ...
        'MarkerFaceColor', [0 0.4470 0.7410], 'MarkerSize', 12)
    axis equal
    axis([-6 -4 -6 -4])
    ax.YTick = [-6:0.5:-4];
    xlabel('log_{10}(S_{true} [m/s])')
    ylabel('log_{10}(S_{est} [m/s])')
    ax.FontSize = 30;
    ax.YAxisLocation = 'right';
    text(-5.9, -4.15, ['\rho = ' num2str(round(param_corr(j,2), 2))], 'FontSize', 30, 'FontWeight', 'bold')
    set(gcf, 'Position', [100 100 2025 2025/2.6667])
end

% Figure S4 - Storativity posterior covariance
figure
clf 
t = tiledlayout(1,3);
for j = 1 : num_case
    ax = nexttile;
    p = pcolor(cgrid{1}, cgrid{2}, reshape(sqrt(post_var(num_cells+1:2*num_cells,j)), num_y, num_x));
    p.LineStyle = 'none';
    hold on
    plot(well_locs(:,1), well_locs(:,2), 'ko', 'LineWidth', 2,...
        'MarkerSize', 12, 'MarkerFaceColor', [0.7529 0 0])
    hold off
    axis equal
    axis square
    axis([-95 95 -95 95])
    ax.XTick = [-80:40:80]; ax.YTick = ax.XTick;
    clim([0.05 0.3])
    if j == 1
        xlabel('X (m)')
        ylabel('Y (m)')
    end
    ax.FontSize = 30;
end
t.TileSpacing = 'loose';
t.Padding = 'loose';
c = colorbar;
c.Layout.Tile = 'east';
c.Label.String = '\sigma log_{10}(S [-])';
c.FontSize = 24;
set(gcf, 'Position', [100 100 2025 2025/2.6667])