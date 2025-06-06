% Oscillatory Hydraulic Tomography

% This code validates the 2D numercial model used to explore the impacts of
% aquifer heterogeneity on oscillatory flow signals. The numerical model is
% validated using the fully confined analytical model developed by
% Rasmussen et al. (2003).

% Code developed by Jeremy Patterson
% Created Dec 2018; Updated Jan 2023

%% Clean Environment
close all; clear; clc

%% Specify Directory
addpath(genpath('/.../.../.../'))

%% Forward Model 1 Setup

%Specify domain 
xmin = -100; xmax = -xmin; dx = 2;
ymin = xmin; ymax = xmax; dy = dx;

domain = struct('x',[],'y',[],'z',[]);
domain.x = [xmin : dx : xmax];
domain.y = [ymin : dy : ymax];
domain.z = [0 1];

numx = numel(domain.x) - 1;
numy = numel(domain.y) - 1;
num_cells = numx * numy;

% Well locations
well_locs = [0 0; ...
             0 20; ...
             20 20; ...
            ];
num_wells = size(well_locs,1);

r = sqrt((well_locs(:,1) - well_locs(1,1)).^2 + (well_locs(:,2) - well_locs(1,2)).^2);

% Grid Domain
[coords, cgrid] = plaid_cellcenter_coord(domain);

% Specify boundary types and boundary values (x / y constant head
% boundaries, no flux in z)
bdry_types = [2; 2; 2; 2; 0; 0];
bdry_vals = zeros(6,1);
bdry_L = 1e-5*ones(6,1);
bdrys = struct('types',bdry_types,'vals',bdry_vals,'leaks',bdry_L);

%% Create Test List
P = logspace(1, 3.25, 10);
V = 1e-2;

test_list = [];
for i = 1:1:numel(P)
    for j = 2 : num_wells
        test_list = [...
            test_list; ...
            (2*pi)/P(i) 1 V*pi/P(i) j ...
            ];
    end
end

%% Create the "official" input files needed by the steady periodic model.

%This step should automatically create the input files needed to run all of
%the steady periodic models. 
lnK = -9.2; lnSs = -11.2;

lnK_true = lnK * ones(num_cells,1);
lnSs_true = lnSs * ones(num_cells,1);
params_true = [lnK_true; lnSs_true];

[inputs] = OHT_create_inputs(well_locs,test_list,domain);
y_fxn = @(params) OHT_run_distribKSs(params, domain, bdrys, inputs, 1);

num_omegas = size(inputs,1);
numobs = size(test_list,1);

y_oht = y_fxn(params_true);
A_syn = y_oht(1:numobs);
B_syn = y_oht(numobs+1:2*numobs);

amp_oht = sqrt(A_syn.^2 + B_syn.^2);
phase_oht = atan2(-B_syn, A_syn);

%% Analytical Solutions
T = exp(lnK);
S = exp(lnSs);
D = T / S;

omega = test_list(:,1);
rad = r(test_list(:,4));
Q_max = test_list(:,3);
for k = 1 : numel(test_list(:,1))
        arg = sqrt(((1i .* rad(k).^2 .* omega(k)) ./ D));
        
        phasor_conf(k,:) = Q_max(k) / (2 * pi * T) * besselk(0, arg);        
        phase_conf(k,:) = -angle(phasor_conf(k));
        amp_conf(k,:) = abs(phasor_conf(k));
end

%% Figures
rad_list = unique(rad);
figure
clf

subplot(1,2,1)
ax = gca;
hold on
plot(P, amp_conf(1:2:end-1), 'k-', 'LineWidth', 2)
plot(P, amp_conf(2:2:end), 'k-', 'LineWidth', 2)
plot(P, amp_oht(1:2:end-1), 'k^', 'MarkerFaceColor', 'k')
plot(P, amp_oht(2:2:end), 'rv', 'MarkerFaceColor', 'r')

xlabel('Period (s)')
ylabel('Amplitude (m)')
ax.FontSize = 18;
ax.XScale = 'log';

subplot(1,2,2)
ax = gca;
hold on
plot(P, phase_conf(1:2:end-1), 'k-', 'LineWidth', 2)
plot(P, phase_conf(2:2:end), 'k-', 'LineWidth', 2)
plot(P, phase_oht(1:2:end-1), 'k^', 'MarkerFaceColor', 'k')
plot(P, phase_oht(2:2:end), 'rv', 'MarkerFaceColor', 'r')
xlabel('Period (s)')
ylabel('Phase Angle (rad)')
ax.FontSize = 18;
ax.XScale = 'log';

set(gcf, 'Position', [100 100 1900 700])