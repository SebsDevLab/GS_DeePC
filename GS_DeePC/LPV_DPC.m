%% File Name: LPV_DPC.m
% Author: Sebastian Zieglmeier
% Date last updated: 30.11.2025
% Description: LPV-DPC, taken from [2], slightly modified,
% Header with originial authors: see below
%
% Sources:
% [1] - Sebastian Zieglmeier, et al., "Gain-Scheduling Data-Enabled 
%       Predictive Control for Nonlinear Systems with Linearized Operating 
%       Regions", 
%       https://doi.org/10.48550/arXiv.2512.02797
% [2] - Chris Verhoek, et al., "Data-Driven Predictive Control for Linear 
%       Parameter-Varying Systems", 
%       https://doi.org/10.1016/j.ifacol.2021.08.588




%% ===================================================================== %%
%  Unbalanced disc example used in the 2021 LPVS paper:
%    C. Verhoek, H. S. Abbas, R. Tóth, and S. Haesaert. "Data-Driven
%    Predictive Control for Linear Parameter-Varying Systems." In Proc. of
%    the 4th IFAC Workshop on Linear Parameter-Varying Systems (LPVS), 2021.
%
%  The code compares a DPC and MPC on a nonlinear system, which is embedded
%  as a self-scheduled LPV system.
%
%  Authors: H. S. Abbas and C. Verhoek
%
%  Released under the BSD 3-Clause License. 
%  Copyright (c) 2022, Eindhoven, The Netherlands
%
%
%  Final version -- Results used in published paper.
%  Tested in ML2020b, ML2021b on both mac and windows machine
%% ===================================================================== %%
clearvars; clc; %close all
rng(350)
plot_results = true;
%% This script compares data-driven lpv predictive controller
% (using certain Hankel matrix) and classic lpv mpc based on an lpvio model
% of an unbalanced disc system


%%  Setting for the predictive controller

nd = 200;
nsim = 50;
ts = 75e-3;

% tuning parameters
Np = 5; % prediction horizon
Q = 1;
R = 0.5e-1;

% constraints

ymax =  pi;
ymin = -pi;
umax =  0.25;
umin = -0.25;

% optimiztion and simulation time
tsim = 20+Np*ts;
tk = (0:ts:tsim)';
ksim = floor(tsim/ts);
tknd = (0:ts:nd)';

% constructing reference signal
r(tk<=4) = pi/2;
r(tk>4 & tk<=8)   = -pi/2;
r(tk>8 & tk<=12)  = pi/2;
r(tk>12 & tk<=16) = 3/4*pi;
r(tk>16 & tk<=20) = -1/4*pi;
r(tk>20) = -1/4*pi;
r = r';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generating LPV data

BAND = [0  0.75];  % w.r.t nyquist freq.
LEVELS = 0.25*[-1 1]; % 1.85 can be change to chnage the power of the signal
SINEDATA = [10, 10, 1];
[u, freq] = idinput([nd,1,1],'sine',BAND, LEVELS, SINEDATA);
x1c = zeros(nd, 1);
x2c = zeros(nd, 1);
x0 = [0 0]';
x0k = x0;
% CT data
for k=1:nd
    [~, x] = ode45(@(t, x) unbaldisc(t, x, u(k)), [0 ts], x0k);
    x0k    = x(end,:)';
    x1c(k+1) = x0k(1);
    x2c(k+1) = x0k(2);
end
y=x1c(1:nd);
p=sin(y)./y;
p(isnan(p))=1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% data-driven lpv predictive control
% Arranging data

nl = 2; % the lag of the system
n  = 2; % minimal order of the system

% the number of data points of u to compute persistently exciting Hankel matrix
%
nh  = nl+Np;
nu  = size(u,2);
ny  = size(y,2);
np  = size(p,2);
pnu = np*nu;
pny = np*ny;
ns  = pnu+nu+pny+ny;
Nd  = (ns+1)*nh-1;
i0  = 1;
id  = Nd+nl+1+i0;



[Hnlu, Hnlup, Hnly, Hnlyp, HLu, HLup, HLy, HLyp] = hanks2(u(i0:id-nl-1,:), y(i0:id-nl-1,:), p(i0:id-nl-1,:), n, nl, Np);

%% ------- Design and implmentation of the proposed lpvdpc

% fixed matrices
zupini = zeros(np*nl, 1);
zypini = zeros(np*nl, 1);
zuph   = zeros(np*Np, 1);
zyph   = zeros(np*Np, 1);

% Settings of the solvar
options = sdpsettings('solver', 'mosek', 'verbose', 0);
g       = sdpvar(size(HLu,2), 1);


% inititalizing the lpvdpc
ydpc  = zeros(ksim-Np, 1);
udpc  = zeros(ksim-Np, 1);
pdpc  = zeros(ksim-Np, 1);
x1dpc = zeros(ksim-Np, 1);
x2dpc = zeros(ksim-Np, 1);
Jdpc  = zeros(ksim-Np, 1);

y0   = pi/2;
x0   = [y0 0]';
x0k  = x0';

uini = zeros(nu*nl,1);
yini = [zeros(nu*(nl-1),1); y0];
pini = sin(yini)./yini;
pini(isnan(pini))=1;

ydpc(1)  = y0;
pdpc(1)  = pini(2);
x1dpc(1) = x0k(1);
x2dpc(1) = x0k(2);

disp('LPV DPC steps:')
for k = 2:ksim-Np
    
    if mod(k,10)==0; fprintf('%.f\n',k); end
    
    [coluini, colyini, dpini] = coldiag(uini, yini, pini);
    
    ph = pini(2)*ones(Np,1);
    
    %  cost function and constraints
    uh = HLu*g;
    yh = HLy*g;
    [coluh, ~, dph] = coldiag(uh, [], ph);
    
    X = [Hnlu
        Hnlup - dpini*Hnlu
        Hnly
        Hnlyp - dpini*Hnly
        HLu
        HLup - dph*HLu
        HLyp - dph*HLy
        ];
    h =[coluini; zupini; colyini; zypini; coluh; zuph; zyph];
    Constraints = X*g == h;
    
    Objective   = 0;
    for ik = 0:Np-1
        Objective   = Objective + (yh(ik+1)-r(k+ik))'*Q*(yh(ik+1)-r(k+ik)) + uh(ik+1)'*R*uh(ik+1);
        Constraints = [Constraints,  umin <= uh(ik+1) <= umax]; %#ok
        Constraints = [Constraints,  ymin <= yh(ik+1) <= ymax]; %#ok
    end
    
    % QP optimization
    diagnostics = optimize(Constraints, Objective ,options);
    opttime = diagnostics.solvertime;
    if diagnostics.problem ~= 0
        % error('infeasible initial condition check the domain of attraction');
    end
    Jdpc(k) = double(Objective);
    gk      = double(g);
    
    % computing control input
    udpc(k) = double(uh(1));
    
    % implementation onto the plant
    [~, x] = ode45(@(t, x) unbaldisc(t, x, udpc(k)), [0 ts], x0k);
    x0k    = x(end,:)';
    x1dpc(k) = x0k(1);
    x2dpc(k) = x0k(2);
    ydpc(k)  = x0k(1);
    
    % initializing the lpvdpc for the next iteration
    uini  = [uini(2:end); udpc(k)];
    yini  = [yini(2:end); ydpc(k)];
    pini  = sin(yini)./yini;
    pini(isnan(pini))=1;
    pdpc(k) = pini(2);
    
end
RMSE_1 = rmse(ydpc(1:ksim-Np), r(1:ksim-Np));
disp(RMSE_1)

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
    subplot(1,3,1); hold on; grid on; box on;set(gca,'fontsize',fontsz_tick)
    plot(tknd(1:nd),y(1:nd),'Color',dgreen, 'LineWidth',linewdt)
    ylabel('$y(k)$ [rad]','Interpreter','latex','fontsize',fontsz_label)
    xlabel('Time [s]','Interpreter','latex','fontsize',fontsz_label)
    dify = 0.05*(max(y(1:nd))-min(y(1:nd)));
    axis([0 tknd(nd) min(y(1:nd))-dify, dify+max(y(1:nd))]);
        
    subplot(1,3,2); hold on; grid on; box on;set(gca,'fontsize',fontsz_tick)
    stairs(tknd(1:nd),u(1:nd),'Color',dgreen, 'LineWidth',linewdt)
    ylabel('$u(k)$ [rad]','Interpreter','latex','fontsize',fontsz_label)
    xlabel('Time [s]','Interpreter','latex','fontsize',fontsz_label)
    axis([0 tknd(nd) -0.28 0.28]);
    
    subplot(1,3,3); hold on; grid on; box on;set(gca,'fontsize',fontsz_tick)
    plot(tknd(1:nd),p(1:nd),'Color',dgreen, 'LineWidth',linewdt)
    ylabel('$p(k)$ [m/rad] ','Interpreter','latex','fontsize',fontsz_label)
    xlabel('Time [s]','Interpreter','latex','fontsize',fontsz_label)
    axis([0 tknd(Nd) 0.89 1.11]); yticklabels({'0.9','','1','','1.1'})

    f2 = figure(2);clf; f2.Position(3:4) = [800 375];    
    subplot(2,2,1);hold on; grid on; box on;set(gca,'fontsize',fontsz_tick)
    plot(tk(1:k),ones(k,1),'Color',gray, 'LineWidth',linewdt, 'HandleVisibility','off')
    plot(tk(1:k),-ones(k,1),'Color',gray, 'LineWidth',linewdt, 'HandleVisibility','off')
    stairs(tk(1:k),r(1:k),'Color',yellow, 'LineWidth',linewdt, 'DisplayName', '$r(k)$')
    plot(tk(1:k),ydpc(1:k),'Color',dgreen, 'LineWidth',linewdt, 'DisplayName', '$y_{\mathrm{DPC}}(k)$')
    ylabel('$y(k)$ [rad]','Interpreter','latex','fontsize',fontsz_label)
    axis([0 tk(k) -1.2 1.2]);
    f2ax1 = gca;
    
    subplot(2,2,2);hold on; grid on; box on;set(gca,'fontsize',fontsz_tick)
    plot(tk(1:k),0.25*ones(k,1),'Color',gray, 'LineWidth',linewdt, 'HandleVisibility','off')
    plot(tk(1:k),-0.25*ones(k,1),'Color',gray, 'LineWidth',linewdt, 'HandleVisibility','off')
    stairs(tk(1:k),udpc(1:k),'Color',dgreen, 'LineWidth',linewdt,'DisplayName','$u_{\mathrm{DPC}}(k)$');
    ylabel('$u(k)$ [V]','Interpreter','latex','fontsize',fontsz_label)
    axis([0 tk(k) -0.3 0.3]);
    f2ax2 = gca;
        
    subplot(2,2,3);hold on; grid on; box on;set(gca,'fontsize',fontsz_tick)
    plot(tk(1:k),ones(k,1),'Color',gray, 'LineWidth',linewdt, 'HandleVisibility','off')
    plot(tk(1:k),pdpc(1:k),'Color',dgreen, 'LineWidth',linewdt,'DisplayName','$p_{\mathrm{DPC}}(k)$');
    ylabel('$p(k)$ [m/rad]','Interpreter','latex','fontsize',fontsz_label)
    axis([0 tk(k) 0.99*min(pdpc(1:k)) 1.01*max(pdpc(1:k))]);
    xlabel({'Time [s]'},'Interpreter','latex','fontsize',fontsz_label)
    f2ax3 = gca;
    
    subplot(2,2,4);hold on; grid on; box on;set(gca,'fontsize',fontsz_tick)
    plot(tk(1:k),x1dpc(1:k),'Color',dgreen, 'LineWidth',linewdt, 'DisplayName', '$x_{1,\mathrm{DPC}}(k)$');
    plot(tk(1:k),x2dpc(1:k),'Color',blue, 'LineWidth',linewdt, 'DisplayName', '$x_{2,\mathrm{DPC}}(k)$');
    ylabel('$x_1(k)$, $x_2(k)$','Interpreter','latex','fontsize',fontsz_label)
    xlabel({'Time [s]'},'Interpreter','latex','fontsize',fontsz_label)
    axis([0 tk(k) -2.5 2.5]);
    f2ax4 = gca;
    
end

%% LPV-MPC control of the unbalanced disc
clear yalmip Objective Constraints J
% settings of DT model 
A0 = [1 ts; 0 (1-ts/0.5971)];
A1 = [0 0; -0.07*9.8*0.42e-3*ts/2.2e-4 0];
B0 = [0; 15.3145*ts/0.5971];

control = sdpvar(Np,1);

% inititalizations
Jmpc = zeros(ksim-Np, 1);
umpc = zeros(ksim-Np, 1);
ympc = zeros(ksim-Np, 1);
pmpc = zeros(ksim-Np, 1);
x1mpc = zeros(ksim-Np, 1);
x2mpc = zeros(ksim-Np, 1);


%%  optimization and simulation
x0k = x0;
y0k = x0(1);
p0k  = sin(y0k)/y0k;
p0k(isnan(p0k)) = 1;

x1mpc(1) = x0k(1);
x2mpc(1) = x0k(2);
ympc(1)  = y0k;
pmpc(1)  = p0k;
umpc(1)  = 0;

disp('Classic LPV MPC steps:')
for k = 2:ksim-Np
    
    if mod(k,10)==0; fprintf('%.f\n',k); end
    %  cost function and constraints
    Objective   = (y0k-r(k))'*Q*(y0k-r(k)) + control(1)'*R*control(1);
    Constraints = umin <= control(1) <= umax; %#ok
    Ap          = A0 + A1*p0k;
    xik         = Ap*x0k + B0*control(1);
    yik         = xik(1);
    for ik = 1:Np-1
        Objective   = Objective   + (yik-r(k+ik))'*Q*(yik-r(k+ik))  + control(ik+1)'*R*control(ik+1);
        Constraints = [Constraints,  umin <= control(ik+1) <= umax]; %#ok 
        Constraints = [Constraints,  ymin <= yik <= ymax]; %#ok
        xik         = Ap*xik + B0*control(ik+1);
        yik         = xik(1);
    end
    
    
    %% QP optimization
    diagnostics = optimize(Constraints, Objective ,options);
    
    if diagnostics.problem ~= 0
        error('infeasible initial condition check the domain of attraction');
    end
    J(k) = double(Objective); %#ok
    Uk   = double(control);
    
    
    %% updates
    % control input update
    umpc(k)     = Uk(1);
    
    % implementation onto the plant
    [~, x] = ode45(@(t, x) unbaldisc(t, x, umpc(k)), [0 ts], x0k);
    x0k    = x(end,:)';
    x1mpc(k) = x0k(1);
    x2mpc(k) = x0k(2);
    y0k      = x0k(1);
    ympc(k) = y0k;
    
    
    % Sch. trajectories update
    p0k = sin(y0k)/y0k;
    p0k(isnan(p0k))=1;
    pmpc(k) = p0k;
end

if plot_results
linewdt = 1.5;
plot(f2ax1,tk(1:k),ympc(1:k),'r--','DisplayName','$y_{\mathrm{MPC}}(k)$', 'LineWidth',linewdt)
stairs(f2ax2,tk(1:k),umpc(1:k),'r--','DisplayName','$u_{\mathrm{MPC}}(k)$', 'LineWidth',linewdt)
plot(f2ax3,tk(1:k),pmpc(1:k),'r--', 'LineWidth',linewdt,'DisplayName','$p_{\mathrm{MPC}}(k)$');
plot(f2ax4,tk(1:k),x1mpc(1:k),'r--', 'LineWidth',linewdt,'DisplayName','$x_{1,\mathrm{MPC}}(k)$');
plot(f2ax4,tk(1:k),x2mpc(1:k),'b--', 'LineWidth',linewdt,'DisplayName','$x_{2,\mathrm{MPC}}(k)$');


legend(f2ax1,'show','fontsize',fontsz_leg);
legend(f2ax2,'show','fontsize',fontsz_leg,'Location','southeast');
legend(f2ax3,'show','fontsize',fontsz_leg,'Location','southeast');
legend(f2ax4,'show','fontsize',fontsz_leg,'Location','best','NumColumns',2);
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
    f2 = figure(3);clf; f2.Position(3:4) = [800 375];    
    subplot(4,1,1);hold on; grid on; box on;set(gca,'fontsize',fontsz_tick)
    plot(tk(1:k),r(1:k),'Color',yellow, 'LineWidth',linewdt, 'DisplayName', '$r(k)$')
    plot(tk(1:k),ydpc(1:k),'Color',dgreen, 'LineWidth',linewdt, 'DisplayName', '$y_{\mathrm{DPC}}(k)$')
    ylabel('$y(k)$ [rad]','Interpreter','latex','fontsize',fontsz_label)
    axis([0 tk(k) -3.5 3.5]);
    f2ax1 = gca;
    
        
    % subplot(4,1,4);hold on; grid on; box on;set(gca,'fontsize',fontsz_tick)
    % plot(tk(1:k), H_idx(1:k), 'LineWidth',linewdt, 'DisplayName','Hankelmatrix Index')
    % xlabel({'Time [s]'},'Interpreter','latex','fontsize',fontsz_label)
    % ylabel('Active Composite Region','Interpreter','latex','fontsize',fontsz_label)
    % f2ax3 = gca;
    
    % subplot(2,2,4);hold on; grid on; box on;set(gca,'fontsize',fontsz_tick)
    % plot(tk(1:k),x1dpc(1:k),'Color',dgreen, 'LineWidth',linewdt, 'DisplayName', '$x_{1,\mathrm{DPC}}(k)$');
    % plot(tk(1:k),x2dpc(1:k),'Color',blue, 'LineWidth',linewdt, 'DisplayName', '$x_{2,\mathrm{DPC}}(k)$');
    % ylabel('$x_1(k)$, $x_2(k)$','Interpreter','latex','fontsize',fontsz_label)
    % xlabel({'Time [s]'},'Interpreter','latex','fontsize',fontsz_label)
    % axis([0 tk(k) -2.5 2.5]);
    % f2ax4 = gca;

    subplot(4,1,2); hold on; grid on; box on; set(gca,'fontsize',fontsz_tick)
    plot(tk(1:k), rel_err, 'Color', red, 'LineWidth', linewdt, ...
         'DisplayName', 'Relative error');
    ylabel('$e_{\mathrm{rel}}(k)$','Interpreter','latex','fontsize',fontsz_label)
    xlabel('Time [s]','Interpreter','latex','fontsize',fontsz_label)

    subplot(4,1,3); hold on; grid on; box on; set(gca,'fontsize',fontsz_tick)
    plot(tk(1:k), rel_err, 'Color', red, 'LineWidth', linewdt, ...
         'DisplayName', 'Relative error');
    ylabel('$e_{\mathrm{rel}}(k)$','Interpreter','latex','fontsize',fontsz_label)
    xlabel('Time [s]','Interpreter','latex','fontsize',fontsz_label)
    axis([0 tk(k) 0 15]);
end

%% ########################################################################
%  #######                                                          #######
%  #######                     LOCAL FUNCTIONS                      #######
%  #######                                                          #######
%  ########################################################################

function dx = unbaldisc(t,x,u) %#ok
% dc motor with an unbalance disc
%
% HSA
% Lübeck 20-01-2021


x1 = x(1);
x2 = x(2);

x1dot = x2;
x2dot = -1/0.5971*x2 - 0.07*9.8*0.42e-3/2.2e-4*sin(x1) + 15.3145/0.5971*u;

dx  = [x1dot; x2dot];
end

function [Hnlu, Hnlup, Hnly, Hnlyp, HLu, HLup, HLy, HLyp] = hanks2(varargin)
% This function constructs the required Hankel matrices for data-
% driven lpv simulation or predictive control.
% It is based on the new formulation in our LPVS paper
%
% Give the full input/output/schduling u/y/p data as column versotrs
% For signals multivariate signals, they are given as matrices as
% s=[s1 s2 ...]
% n: order of the system
% nl: the lag of the system, which is related to the l-step
%     observabaility matrix Q_l, l that enables us to determine x(0) or
%     that makes Q_l has rank n
% L  %%<<<========= Specify here the number of k-step a head for
%                   simulation or prediction horizon of predictive control
%
% H. S. Abbas
% Lübeck 11.02.2021
%
% NOTE: Stripped down version...


%% Input parsing
u  = varargin{1};
y  = varargin{2};
p  = varargin{3};
%n  = varargin{4};
nl = varargin{5};
L  = varargin{6};
N = size(u,1);

%% constructing io signals with scheduling parameter (only works for ny=nu=1)
up = p.*u;
yp = p.*y;

%% Hankel matrices
Hnlu = hank2(u(1:N-L,:), nl);
Hnlup = hank2(up(1:N-L,:), nl);
Hnly = hank2(y(1:N-L,:), nl);
Hnlyp = hank2(yp(1:N-L,:), nl);
HLu = hank2(u(nl+1:N,:), L);
HLup = hank2(up(nl+1:N,:), L);
HLy = hank2(y(nl+1:N,:), L);
HLyp = hank2(yp(nl+1:N,:), L);
end

function H = hank2(u, L)
% This function constructs the Hankel matrix for u
% u: is a colum vector or matrix consists of colum vectors of signals
% L: No. row blocks in the resulted Hankel matrix
% The update of the function hank where now we use matlab function 'hankel'
%
% HSA
% Lübeck 11.2.2021

u = u';
[m, n] = size(u);
i = hankel(1:n);
i = i(1:L, 1:n-L+1);
H = reshape(u(:,i), m*L, []);
end

function [colu, coly, diagp] = coldiag(u, y, p)
% Constructing columns from u, y input/output, respc. and a block diagonal
% matrix from p
% Multivariate signals should be as s=[s1 s2 ...], si are column vectors
%
% HSA
% Lübeck 14.02.2021


%% setting up dimensions
[N, nu]  = size(u);
if ~isempty(y)
    ny  = size(y,2);
end



%% the column signals correspnding to the initial conditions uini, yini, pini
colu   = reshape(u',[nu*N,1]);
if ~isempty(y)
    coly   = reshape(y',[ny*N,1]);
else
    coly   = [];
end
p  = p';
Cp = mat2cell(p, size(p,1), ones(1,size(p,2)));    %// step 1
diagp = blkdiag(Cp{:});
end

