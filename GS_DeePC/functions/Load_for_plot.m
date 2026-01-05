%% File Name: Load_for_plot.m
% Author: Sebastian Zieglmeier
% Date last updated: 30.11.2025
% Description: Load all result files for comparison plot
%
% Sources:
% [1] - Sebastian Zieglmeier, et al., "Gain-Scheduling Data-Enabled 
%       Predictive Control for Nonlinear Systems with Linearized Operating 
%       Regions", 
%       https://doi.org/10.48550/arXiv.2512.02797
%  
%  This script loads the results of different GS-DeePC variants
%  into separate workspace variables.

clear all; clc;

fprintf("=== Loading result files for comparison ===\n\n");

%% === Load GS-DeePC (n = 16) ======================================
try
    load results_GS_DeePC_n16_m007_lini_1e6_lg_1e3.mat
    R_GS16 = results;
    fprintf("Loaded: results_GS_DeePC_n16_m007_lini_1e6_lg_1e3\n");
catch
    warning("Could not load results_GS_DeePC_n16_m007_lini_1e6_lg_1e3.mat");
end

%% === Load GS-DeePC (n = 8) =======================================
try
    load results_GS_DeePC_n8_m007_lini_1e6_lg_1e3.mat
    R_GS8 = results;
    fprintf("Loaded: results_GS_DeePC_n8_m007_lini_1e6_lg_1e3\n");
catch
    warning("Could not load results_GS_DeePC_n8_m007_lini_1e6_lg_1e3.mat");
end

%% === Load GS-DeePC (n = 4) =======================================
try
    load results_GS_DeePC_n4_m007_lini_1e7_lg_1e3.mat
    R_GS4 = results;
    fprintf("Loaded: results_GS_DeePC_n4_m007_lini_1e7_lg_1e3\n");
catch
    warning("Could not load results_GS_DeePC_n4_m007_lini_1e7_lg_1e3.mat");
end

%% === Load standard DeePC =========================================
try
    load results_DeePC_m007_lini_1e6_lg_1e3.mat
    R_DeePC = results;
    fprintf("Loaded: results_DeePC_m007_lini_1e6_lg_1e3\n");
catch
    warning("Could not load results_DeePC_m007_lini_1e6_lg_1e3.mat");
end

%% === Load LPV-DPC =================================================
try
    load results_LPV_DPC_m007.mat
    R_LPV = results;
    fprintf("Loaded: results_LPV_DPC_m007\n");
catch
    warning("Could not load results_LPV_DPC_m007.mat");
end

%% === Summary =======================================================
fprintf("\n=== Completed loading all available result files ===\n");
fprintf("Variables now available:\n");
fprintf("  R_GS16, R_GS8, R_GS4, R_DeePC, R_LPV\n\n");
fprintf("Next step: run the plotting script.\n");


%% ===============================================================
% HIGH-QUALITY COMPARISON PLOT FOR MULTIPLE DPC APPROACHES
% Creates 3 subplots:
%   (1) Reference tracking
%   (2) Relative error
%   (3) Active region index H_idx
%
% Adjust font sizes, line widths, colors as needed.
% ===============================================================

clear fig;

%% === FIGURE SETTINGS ============================================
fontsz_label = 16;
fontsz_tick  = 14;
linewdt      = 2;

fig = figure(100); clf;
fig.Color = 'w';         % white background
fig.Position(3:4) = [1200 900];   % adjust size manually afterwards if needed

%% === Extract signals shorter of all runs ==========================
% (Ensures same length if slight differences exist)
len = min([ length(R_GS16.tk), length(R_GS8.tk), length(R_GS4.tk), ...
            length(R_DeePC.tk), length(R_LPV.tk) ]);

t  = R_GS16.tk(1:len);
r  = R_GS16.r(1:len);   % reference is identical across approaches

%% Collect outputs
y_GS16 = R_GS16.ydpc(1:len);
y_GS8  = R_GS8.ydpc(1:len);
y_GS4  = R_GS4.ydpc(1:len);
y_DeePC = R_DeePC.ydpc(1:len);
y_LPV   = R_LPV.ydpc(1:len);

%% Relative errors
e_GS16 = R_GS16.rel_err(1:len);
e_GS8  = R_GS8.rel_err(1:len);
e_GS4  = R_GS4.rel_err(1:len);
e_DeePC = R_DeePC.rel_err(1:len);
e_LPV   = R_LPV.rel_err(1:len);

%% Active region indices
h_GS16 = R_GS16.H_idx(1:len);
h_GS8  = R_GS8.H_idx(1:len);
h_GS4  = R_GS4.H_idx(1:len);
% DeePC + LPV have no region switching → set dummy values
h_DeePC = zeros(len,1);
h_LPV   = zeros(len,1);

%% Colors (adjust as desired)
col_ref  = [237,192,1]/255;      % yellow
col1 = [0 75 90]/255;            % GS16 dark green
col2 = [34,139,34]/255;          % GS8 green
col3 = [0, 0, 128]/255;          % GS4 blue
col4 = [181 22 33]/255;          % DeePC dark red
col5 = [0.2 0.2 0.2];            % LPV gray

%% ===============================================================
% SUBPLOT 1 — REFERENCE TRACKING
% ===============================================================

ax1=subplot(3,1,1); hold on; grid on; box on;
set(gca,'fontsize',fontsz_tick)

plot(t, r,      'Color', col_ref, 'LineWidth', linewdt, ...
     'DisplayName', 'Reference');

plot(t, y_GS16, 'Color', col1,    'LineWidth', linewdt, ...
     'DisplayName', 'GS-DeePC M=16');

plot(t, y_GS8,  'Color', col2,    'LineWidth', linewdt, ...
     'DisplayName', 'GS-DeePC M=8');

plot(t, y_GS4,  'Color', col3,    'LineWidth', linewdt, ...
     'DisplayName', 'GS-DeePC M=4');

plot(t, y_DeePC,'Color', col4,    'LineWidth', linewdt, ...
     'DisplayName', 'DeePC');

plot(t, y_LPV,  'Color', col5,    'LineWidth', linewdt, ...
     'DisplayName', 'LPV-DPC');

ylabel('$y(k)$ [rad]','Interpreter','latex','fontsize',fontsz_label)
lgd = legend('Location','best')   % dummy legend to modify later
lgd.FontSize = 12;
lgd.NumColumns = 3;
title('')                   % no title

%% ===============================================================
% ZOOM-IN INSET FOR SUBPLOT 1
% ===============================================================

% --- Coordinates of zoom window ---
x_zoom = [9 12];
y_zoom = [1.5 1.7];

% --- Create inset axes ---
inset_width  = 0.18;   % adjust inset width
inset_height = 0.10;   % adjust inset height
inset_xpos   = 0.55;   % relative position (0–1)
inset_ypos   = 0.65;   % relative position (0–1)

ax_inset = axes('Position', [inset_xpos inset_ypos inset_width inset_height]);
box on; hold on; grid on;

% --- Plot same signals inside inset ---
plot(t, r,      'Color', col_ref, 'LineWidth', 1.2);
plot(t, y_GS16, 'Color', col1,    'LineWidth', 1.2);
plot(t, y_GS8,  'Color', col2,    'LineWidth', 1.2);
plot(t, y_GS4,  'Color', col3,    'LineWidth', 1.2);
plot(t, y_DeePC,'Color', col4,    'LineWidth', 1.2);
plot(t, y_LPV,  'Color', col5,    'LineWidth', 1.2);

% --- Set zoom limits ---
xlim(x_zoom);
ylim(y_zoom);

set(ax_inset, 'FontSize', fontsz_tick - 2);

%% ===============================================================
% SUBPLOT 2 — RELATIVE ERROR
% ===============================================================

ax2 = subplot(3,1,2); hold on; grid on; box on;
set(gca,'fontsize',fontsz_tick)

plot(t, e_GS16*100, 'Color', col1, 'LineWidth', linewdt, 'DisplayName','GS16 dummy');
plot(t, e_GS8*100,  'Color', col2, 'LineWidth', linewdt, 'DisplayName','GS8 dummy');
plot(t, e_GS4*100,  'Color', col3, 'LineWidth', linewdt, 'DisplayName','GS4 dummy');
plot(t, e_DeePC*100,'Color', col4, 'LineWidth', linewdt, 'DisplayName','DeePC dummy');
plot(t, e_LPV*100,  'Color', col5, 'LineWidth', linewdt, 'DisplayName','LPV dummy');

ylabel('$e_{\mathrm{rel,ss}}(k)$ [\%]','Interpreter','latex','fontsize',fontsz_label)

% legend('Location','best')   % dummy legend
title('')

%% ===============================================================
% SUBPLOT 3 — ACTIVE REGION INDEX
% ===============================================================

ax3 = subplot(3,1,3); hold on; grid on; box on;
set(gca,'fontsize',fontsz_tick)

plot(t, h_GS16, 'Color', col1, 'LineWidth', linewdt, 'DisplayName','GS16 dummy');
plot(t, h_GS8,  'Color', col2, 'LineWidth', linewdt, 'DisplayName','GS8 dummy');
plot(t, h_GS4,  'Color', col3, 'LineWidth', linewdt, 'DisplayName','GS4 dummy');
plot(t, h_DeePC,'Color', col4, 'LineWidth', linewdt, 'DisplayName','DeePC dummy');
plot(t, h_LPV,  'Color', col5, 'LineWidth', linewdt, 'DisplayName','LPV dummy');

ylabel('$C_{\mathrm{i}}$ [-]','Interpreter','latex','fontsize',fontsz_label)
xlabel('Time [s]','Interpreter','latex','fontsize',fontsz_label)
% legend('Location','best')
title('')

% Desired geometry (normalized figure units):
% Heights in ratio 1 : 1 : 0.5 → e.g. 0.32, 0.32, 0.16
h1 = 0.3;
h2 = 0.2;
h3 = 0.2;
gap = 0.05;      % small vertical gap between axes
bottom3 = 0.08;  % bottom margin

y3 = bottom3;
y2 = y3 + h3 + gap;
y1 = y2 + h2 + gap;

left  = 0.12;
width = 0.82;

set(ax1, 'Position', [left y1 width h1]);
set(ax2, 'Position', [left y2 width h2]);
set(ax3, 'Position', [left y3 width h3]);

exportgraphics(fig, 'comparison_plot.png', 'Resolution', 600);

