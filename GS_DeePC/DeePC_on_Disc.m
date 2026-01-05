%% File Name: DeePC_on_Disc.m
% Author: Sebastian Zieglmeier
% Date last updated: 30.11.2025
% Description: Standard DeePC on nonlinear unbalanced disc system to find
% limitations of standard DeePC.
%
% Sources:
% [1] - Sebastian Zieglmeier, et al., "Gain-Scheduling Data-Enabled 
%       Predictive Control for Nonlinear Systems with Linearized Operating 
%       Regions", 
%       https://doi.org/10.48550/arXiv.2512.02797
% [2] - Chris Verhoek, et al., "Data-Driven Predictive Control for Linear 
%       Parameter-Varying Systems", 
%       https://doi.org/10.1016/j.ifacol.2021.08.588
%
% Functions:
% DeePC_fast_GS.m
% make_hankel.m, make_page.m


clearvars; close all; clc

%%  Setting for the predictive controller
plot_results = 1;
nd = 400;             % number data points
Tsim = 20;              % Simulation Time
ts = 75e-3;             % Sampling Time

% tuning parameters
Q = 10;
R = 0.5e-1;
T_ini = 2;              % initialization horizon
T_fut = 5;              % prediction horizon
lambda_ini = 1e5;
lambda_g = 1e2;

% constraint and system parameters
sys.constraints.u_min = -0.25;
sys.constraints.u_max = 0.25;
sys.constraints.y_min = -pi;
sys.constraints.y_max = pi;
sys.nu = 1; 
sys.ny = 1; 

% optimiztion and simulation time
tsim = Tsim+T_fut*ts;
tk = (0:ts:tsim)';
tknd = (0:ts:nd)';
ksim = floor(tsim/ts);


% constructing reference signal
r(tk<=4) = pi/2;
r(tk>4 & tk<=8)   = -pi/2;
r(tk>8 & tk<=12)  = pi/2;
r(tk>12 & tk<=16) = 3/4*pi;
r(tk>16 & tk<=20) = -1/4*pi;
r(tk>20) = -1/4*pi;
r = r';

%% Generating sufficient data spread over all regions
u = zeros(nd, 1);
x1c = zeros(nd, 1);
x2c = zeros(nd, 1);
x0 = [0.25 0]';
x0k = x0;
upper_limit = 0.9*sys.constraints.y_max;
lower_limit = 0.9*sys.constraints.y_min;
% higher and lower mean to bang up or down when reaching the overall system
% limits:
bang_down = 0;
bang_up = 0.05;
u_ripple = 0.25*0.2*idinput([nd,1,1], 'RBS');
u_ripple2 = 0.25*0.2*idinput([nd,1,1], 'RBS');
u_ripple3 = 0.25*0.2*idinput([nd,1,1], 'RBS');
u_ripple4 = 0.25*0.2*idinput([nd,1,1], 'RBS');
% Collect data
for k=1:nd
    if x0k(1) >= upper_limit
        bang_down = 0.05;
        bang_up = 0;
    elseif x0k(1) <= lower_limit
        bang_up = 0.05;
        bang_down = 0;
    else
        % do nothing
    end
    u(k) = u_ripple(k) + u_ripple2(k) + u_ripple3(k) + u_ripple4(k) + (bang_up-bang_down);
    u(k) = max(u(k), sys.constraints.u_min);
    u(k) = min(u(k), sys.constraints.u_max);
    [~, x] = ode45(@(t, x) unbaldisc(t, x, u(k)), [0 ts], x0k);
    x0k    = x(end,:)';
    x1c(k+1) = x0k(1);
    x2c(k+1) = x0k(2);
end
y=x1c(1:nd);
u_data = u;
y_data = y;

%% Initilaize Controller
T_d = length(u_data);
L = T_ini + T_fut; % Lag for right size of Hankel matrizes
num_hankel_cols = T_d - L + 1;
H_u = make_hankel_MIMO(u_data',num_hankel_cols, L);
H_y = make_hankel_MIMO(y_data',num_hankel_cols, L);

control = DeePC_fast(T_d, T_ini, T_fut, R, Q, (lambda_ini), lambda_g, H_u, H_y, sys);

%% Initialize Simulation
% start conditions
ydpc  = zeros(ksim-T_fut, 1);
udpc  = zeros(ksim-T_fut, 1);
x1dpc = zeros(ksim-T_fut, 1);
x2dpc = zeros(ksim-T_fut, 1);
H_idx  = zeros(ksim-T_fut, 1);

y0 = r(1);
x0k   = [y0 0];

uini = zeros(sys.nu*T_ini,1);
yini = [zeros(sys.ny*(T_ini-1),1); y0];
ydpc(1)  = y0;
x1dpc(1) = x0k(1);
x2dpc(1) = x0k(2);

for k = 2:ksim-T_fut
    % reference for the prediction horizon
    y_reference = r(k:k+T_fut-1, :);
    %  Call the DeePC-GS controller step with selected Hankel matrices 
    [u_fut, y_fut, g_fut, sigma_fut, opt_error_flag] = control.step(uini, yini, y_reference, H_u, H_y); 
    if opt_error_flag
        % throw warning and use fallback: u=0
        warning('k=%d: optimization returned error flag. Using u=0 fallback.', k);
        udpc(k) = 0;
    else
        udpc(k) = u_fut(1);
    end

    %  apply to the system (simulate one sampling period) 
    [~, x] = ode45(@(t, x) unbaldisc(t, x, udpc(k)), [0 ts], x0k);
    x0k    = x(end,:)';
    x1dpc(k) = x0k(1);
    x2dpc(k) = x0k(2);
    ydpc(k)  = x0k(1);

    %  update Hankel/past stacks for next iteration 
    uini  = [uini(2:end); udpc(k)];
    yini  = [yini(2:end); ydpc(k)];
end


function dx = unbaldisc(t,x,u) 
% dc motor with an unbalance disc from [2]
x1 = x(1);
x2 = x(2);

x1dot = x2;
x2dot = -1/0.5971*x2 - 0.07*9.8*0.42e-3/2.2e-4*sin(x1) + 15.3145/0.5971*u;

dx  = [x1dot; x2dot];
end



% Numerical results: 
%% Steady-state intervals (seconds)
intervals = [ 2  3;
              6  7;
              10 11;
             14 15;
             18 19 ];

%% Compute pointwise relative error
eps_rel = 1e-1;                % avoid division by zero
rel_err_full  = abs(ydpc(1:k) - r(1:k)) ./ max(eps_rel, abs(r(1:k)));
abs_err_full = abs(ydpc(1:k) - r(1:k));
%% Build mask for steady-state intervals only
mask = false(size(tk(1:k)));
for i = 1:size(intervals,1)
    t1 = intervals(i,1);
    t2 = intervals(i,2);
    mask = mask | (tk(1:k) >= t1 & tk(1:k) <= t2);
end

%% Relative error signal (zero outside intervals)
rel_err = zeros(size(rel_err_full));
rel_err(mask) = rel_err_full(mask);

%% Compute RMSE only on steady-state intervals
rmse_ss = sqrt(mean(abs_err_full(mask).^2));
rmse_trans = sqrt(mean(abs_err_full(~mask).^2));

%% Display
fprintf("RMSE of steady-state error = %.6f\n", rmse_ss);
fprintf("RMSE of transient error = %.6f\n", rmse_trans);


%% plots
% plot settings
fontsz_label = 16;
fontsz_tick = 12;
fontsz_leg = 14;
linewdt = 2;
% colors
dgreen = [0 75 90]/255;
green   = [0, 0, 128]/255;
blue  = [34,139,34]/255;
dred    = [181 22 33]/255;
red   = [255 22 33]/255;
yellow = [237,192,1]/255;
gray   = [1 1 1]*0.75;
black   = [1 1 1]*0;

%plot
if plot_results

    f1 = figure(1);clf; f1.Position(3:4) = [800 150];
    subplot(1,2,1); hold on; grid on; box on;set(gca,'fontsize',fontsz_tick)
    plot(tknd(1:nd),y(1:nd),'Color',dgreen, 'LineWidth',linewdt)
    ylabel('$y(k)$ [rad]','Interpreter','latex','fontsize',fontsz_label)
    xlabel('Time [s]','Interpreter','latex','fontsize',fontsz_label)
    dify = 0.05*(max(y(1:nd))-min(y(1:nd)));
    axis([0 tknd(nd) min(y(1:nd))-dify, dify+max(y(1:nd))]);
        
    subplot(1,2,2); hold on; grid on; box on;set(gca,'fontsize',fontsz_tick)
    stairs(tknd(1:nd),u(1:nd),'Color',dgreen, 'LineWidth',linewdt)
    ylabel('$u(k)$ [rad]','Interpreter','latex','fontsize',fontsz_label)
    xlabel('Time [s]','Interpreter','latex','fontsize',fontsz_label)
    axis([0 tknd(nd) -0.28 0.28]);
    


    f2 = figure(2);clf; f2.Position(3:4) = [800 375];    
    subplot(2,2,1);hold on; grid on; box on;set(gca,'fontsize',fontsz_tick)
    plot(tk(1:k),r(1:k),'Color',yellow, 'LineWidth',linewdt, 'DisplayName', '$r(k)$')
    plot(tk(1:k),ydpc(1:k),'Color',dgreen, 'LineWidth',linewdt, 'DisplayName', '$y_{\mathrm{DPC}}(k)$')
    ylabel('$y(k)$ [rad]','Interpreter','latex','fontsize',fontsz_label)
    axis([0 tk(k) -3.5 3.5]);
    f2ax1 = gca;
    
    subplot(2,2,2);hold on; grid on; box on;set(gca,'fontsize',fontsz_tick)
    plot(tk(1:k),0.25*ones(k,1),'Color',gray, 'LineWidth',linewdt, 'HandleVisibility','off')
    plot(tk(1:k),-0.25*ones(k,1),'Color',gray, 'LineWidth',linewdt, 'HandleVisibility','off')
    stairs(tk(1:k),udpc(1:k),'Color',dgreen, 'LineWidth',linewdt,'DisplayName','$u_{\mathrm{DPC}}(k)$');
    ylabel('$u(k)$ [V]','Interpreter','latex','fontsize',fontsz_label)
    axis([0 tk(k) -0.3 0.3]);
    f2ax2 = gca;
        
    subplot(2,2,4);hold on; grid on; box on;set(gca,'fontsize',fontsz_tick)
    plot(tk(1:k), H_idx(1:k), 'LineWidth',linewdt, 'DisplayName','Hankelmatrix Index')
    xlabel({'Time [s]'},'Interpreter','latex','fontsize',fontsz_label)
    f2ax3 = gca;
    

    subplot(2,2,3); hold on; grid on; box on; set(gca,'fontsize',fontsz_tick)
    plot(tk(1:k), rel_err, 'Color', red, 'LineWidth', linewdt, ...
         'DisplayName', 'Relative error');
    ylabel('$e_{\mathrm{rel}}(k)$','Interpreter','latex','fontsize',fontsz_label)
    xlabel('Time [s]','Interpreter','latex','fontsize',fontsz_label)

    f2 = figure(3);clf; f2.Position(3:4) = [800 375];    
    subplot(4,1,1);hold on; grid on; box on;set(gca,'fontsize',fontsz_tick)
    plot(tk(1:k),r(1:k),'Color',yellow, 'LineWidth',linewdt, 'DisplayName', '$r(k)$')
    plot(tk(1:k),ydpc(1:k),'Color',dgreen, 'LineWidth',linewdt, 'DisplayName', '$y_{\mathrm{DPC}}(k)$')
    ylabel('$y(k)$ [rad]','Interpreter','latex','fontsize',fontsz_label)
    axis([0 tk(k) -3.5 3.5]);
    f2ax1 = gca;
    
    subplot(4,1,4);hold on; grid on; box on;set(gca,'fontsize',fontsz_tick)
    plot(tk(1:k), H_idx(1:k), 'LineWidth',linewdt, 'DisplayName','Hankelmatrix Index')
    xlabel({'Time [s]'},'Interpreter','latex','fontsize',fontsz_label)
    ylabel('Active Composite Region','Interpreter','latex','fontsize',fontsz_label)
    f2ax3 = gca;
    
    subplot(4,1,2); hold on; grid on; box on; set(gca,'fontsize',fontsz_tick)
    plot(tk(1:k), rel_err*100, 'Color', red, 'LineWidth', linewdt, ...
         'DisplayName', 'Relative error');
    ylabel('$e_{\mathrm{rel}}(k)$','Interpreter','latex','fontsize',fontsz_label)
    xlabel('Time [s]','Interpreter','latex','fontsize',fontsz_label)

    subplot(4,1,3); hold on; grid on; box on; set(gca,'fontsize',fontsz_tick)
    plot(tk(1:k), rel_err*100, 'Color', red, 'LineWidth', linewdt, ...
         'DisplayName', 'Relative error');
    ylabel('$e_{\mathrm{rel}}(k)$','Interpreter','latex','fontsize',fontsz_label)
    xlabel('Time [s]','Interpreter','latex','fontsize',fontsz_label)
    axis([0 tk(k) 0 15]);
end
