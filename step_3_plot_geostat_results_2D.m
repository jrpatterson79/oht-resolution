% Plot geostatistical imaging results

%% Clean Environment
close all; clear; clc

%% Specify Directory
load_dir = '/.../.../'; % Directory with inversion output .mat files 
print_dir = '/.../.../';
addpath(genpath('/.../.../'))
%% Load File
% Uncomment these lines to plot geostatistical inversions     
file_name = {'geostat_1_freq';
             'geostat_2_freq';
             'geostat_4_freq';
             };
num_case = numel(file_name);

for w = 1 : num_case
    reg_idx = 9;
    load([load_dir file_name{w}])
    disp(file_name{w})
    disp(['K_var = ' num2str(lnK_var_est(reg_idx)) '; Ss_var = ' num2str(lnSs_var_est(reg_idx))])

    param_corr(w,:) = corr(params_true(coords_range), params_best_mat(coords_range,reg_idx));
    params_plot(:,w) = log10(exp(params_best_mat(1:num_cells,reg_idx)));
    K_range = log10(exp((lnK_range)));

    %% Calculate Posterior Covariance (K only) 
    reg_param = lnK_var_est(reg_idx); % K variance

    XK = X(1:num_cells,1);
    H_local = H_tilde_func(params_best_mat(:,reg_idx));
    HK = [H_local(1:num_obs,1:num_cells); H_local(num_obs+1:2*num_obs,1:num_cells)];

    num_drift = size(XK,2);
    QHt = zeros(num_cells, 2*num_obs);
    Q_row_reg = Q_row_corr .* reg_param;
    Qprod_func = @(invec) covar_product_K(Q_row_corr.*lnK_var_est(reg_idx),invec,num_x,num_y);

    for k = 1 : 2*num_obs
        QHt(:,k) = Qprod_func(HK(k,:)');
    end

    HQHt = [HK*QHt + R];
    HX = HK * XK;
    QHtX = [QHt XK];

    covmult_part = [HQHt HX; HX' zeros(num_drift,num_drift)] \ QHtX';
    prior_var = Q_row_reg(1) .* ones(num_cells,1);

    for i = 1 : num_cells
        % Posterior variance matrix - hydraulic conductivity
        post_var(i,w) = log10(exp(prior_var(i) - QHtX(i,:)*covmult_part(:,i)));
    end 
    % Normalized error ((K_mod-K_true)/post_variance)
    norm_err(:,w) = (params_plot(:,w) - log10(exp(params_true(1:num_cells)))) ./ sqrt(post_var(:,w));
end

%% Plot Figures
% Figure 4 - Geostatistical Tomogram with parameter cross-plot
labels = {'1 Frequency'; '2 Frequencies'; '4 Frequencies'};
figure(4)
clf
for j =  1 : num_case
    subplot(3,2,(2*j)-1)
    ax = gca;
    p = pcolor(cgrid{1}, cgrid{2}, reshape(params_plot(:,j), num_y, num_x));
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
    caxis([K_range(1) K_range(2)])
    c = colorbar;
    c.Label.String = 'log_{10}(K_{est} [m/s])';
    c.FontSize = 28;

    subplot(3,2,2*j)
    ax = gca;
    plot([-6 -2], [-6 -2], '-', 'Color', 0.5*[1 1 1], 'LineWidth', 3)
    hold on
    plot(log10(exp(params_true(coords_range))), params_plot(coords_range,j), 'ko', ...
        'MarkerFaceColor', [0 0.4470 0.7410], 'MarkerSize', 12)
    axis equal
    axis([-6 -2 -6 -2])
    ax.YTick = [-6:1:-2];
    xlabel('log_{10}(K_{true} [m/s])')
    ylabel('log_{10}(K_{est} [m/s])')
    ax.FontSize = 30;
    ax.YAxisLocation = 'right';
    text(-5.8, -2.3, ['\rho = ' num2str(round(param_corr(j), 2))], 'FontSize', 30, 'FontWeight', 'bold')
end
set(gcf, 'Position', [100 100 1440 1440])

% Figure 5 - Posterior covariance
figure(5)
clf 
t = tiledlayout(1,3);
for j = 1 : num_case
    ax = nexttile;
    p = pcolor(cgrid{1}, cgrid{2}, reshape(sqrt(post_var(:,j)), num_y, num_x));
    p.LineStyle = 'none';
    hold on
    plot(well_locs(:,1), well_locs(:,2), 'ko', 'LineWidth', 2,...
        'MarkerSize', 12, 'MarkerFaceColor', [0.7529 0 0])
    hold off
    axis equal
    axis square
    axis([-95 95 -95 95])
    ax.XTick = [-80:40:80]; ax.YTick = ax.XTick;
%     caxis([0.1 0.3])
    caxis([0.03 0.3])
    if j == 1
        xlabel('X (m)')
        ylabel('Y (m)')
    end
    ax.FontSize = 30;
    title(labels{j}, 'FontWeight', 'bold', 'FontSize', 30)
end
t.TileSpacing = 'loose';
t.Padding = 'loose';
c = colorbar;
c.Ticks = [0.03 0.09 0.16 0.23 0.3];
c.Layout.Tile = 'east';
c.Label.String = '\sigma log_{10}(K [m/s])';
c.FontSize = 24;
set(gcf, 'Position', [100 100 2025 2025/2.6667])