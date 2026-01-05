%% File Name: GS_DeePC_overlapping_looped.m
% Author: Sebastian Zieglmeier
% Date last updated: 30.11.2025
% Description: Gain Scheduling DeePC with overlapping data-matrices for 
% robustifying DeePC for control of nonlinear systems looped to see
% dependency of overlapping region numbers on the results
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
%
% GS_DeePC_overlapping looped for 50 times with region numbers from 1 to 50
% to measure the dependency of the results on the total number of 
% linearization regions 

clearvars; close all; clc
%%  Setting for the predictive controller
% Result vector intialization:
RMSE_SS = [];
RMSE_trans = [];
number_range = [];

plot_results = 1;
nd = 10000000;          % number data points
Tsim = 20;              % Simulation Time
ts = 75e-3;             % Sampling Time

% DeePC tuning parameters
T_fut = 5;
Q = 100;
R = 0.5e-1;
T_ini = 2;              % initialization horizon
T_fut = 5;              % prediction horizon
lambda_ini = 1e6;
lambda_g = 1e3;

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
x0k = [0.25 0]';
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

%% Start big dependency loop:
for nn=2:1:50
%% Splitting data into regions:
num_ranges = nn;
col_cut = 200;
range_edges = linspace(sys.constraints.y_min, sys.constraints.y_max, num_ranges + 1);
min_segment_length = (T_ini + T_fut);   % minimum segment length to stay inside a region
L = T_ini + T_fut;                      % total Hankel window length
num_hankel_cols = 1;                    % dummy, not needed
num_page_cols = 1;                      % not used in make_hankel_MIMO but can remain for clarity

%% Initialize region structure
regions = struct();
for rr = 1:num_ranges
    regions(rr).range = [range_edges(rr), range_edges(rr+1)];
    regions(rr).segments_y = {};
    regions(rr).segments_u = {};
    regions(rr).H_u = [];  % overall Hankel matrix for u
    regions(rr).H_y = [];  % overall Hankel matrix for y
end

%% Step 1: Segment data by operating range
current_range = NaN;
start_idx = 1;

for k = 1:length(y_data)
    % Determine active range
    in_range = NaN;
    for rr = 1:num_ranges
        lower = range_edges(rr);
        upper = range_edges(rr+1);
        if (rr < num_ranges && y_data(k) >= lower && y_data(k) < upper) || ...
           (rr == num_ranges && y_data(k) >= lower && y_data(k) <= upper)
            in_range = rr;
            break;
        end
    end

    % Continue within same range
    if k > 1 && in_range == current_range
        continue;
    end

    % Handle previous segment
    if k > 1 && ~isnan(current_range)
        end_idx = k - 1;
        segment_len = end_idx - start_idx + 1;
        if segment_len >= min_segment_length
            y_segment = y_data(start_idx:end_idx).';
            u_segment = u_data(start_idx:end_idx).';
            regions(current_range).segments_y{end+1} = y_segment;
            regions(current_range).segments_u{end+1} = u_segment;
        end
    end

    % Update range
    if ~isnan(in_range)
        current_range = in_range;
        start_idx = k;
    else
        current_range = NaN;
    end
end

% Handle final open segment
if ~isnan(current_range)
    segment_len = length(y_data) - start_idx + 1;
    if segment_len >= min_segment_length
        y_segment = y_data(start_idx:end).';
        u_segment = u_data(start_idx:end).';
        regions(current_range).segments_y{end+1} = y_segment;
        regions(current_range).segments_u{end+1} = u_segment;
    end
end

%% Step 2: Build Hankel matrices per segment and combine them
for rr = 1:num_ranges
    all_Hu = [];
    all_Hy = [];
    for s = 1:numel(regions(rr).segments_y)
        u_seg = regions(rr).segments_u{s};
        y_seg = regions(rr).segments_y{s};

        % Ensure row vectors for MIMO consistency
        if iscolumn(u_seg), u_seg = u_seg.'; end
        if iscolumn(y_seg), y_seg = y_seg.'; end

        % Call Page Matrix builder
        % num_page_cols = floor(size(u_seg, 2)/L);
        % H_u = make_page_MIMO(u_seg, num_page_cols, L);
        % H_y = make_page_MIMO(y_seg, num_page_cols, L);
        
        % Call Hankel Matrix builder
        H_u = make_hankel_MIMO(u_seg, num_hankel_cols, L);
        H_y = make_hankel_MIMO(y_seg, num_hankel_cols, L);

        % Concatenate all Hankel/Page matrices
        all_Hu = [all_Hu, H_u];
        all_Hy = [all_Hy, H_y];
    end

    % Store overall Hankel/Page matrices for this region
    regions(rr).H_u = all_Hu;
    regions(rr).H_y = all_Hy;

    fprintf('Region %d: built %d Hankel columns (u: %dx%d, y: %dx%d)\n', ...
        rr, size(all_Hu,2), size(all_Hu,1), size(all_Hu,2), size(all_Hy,1), size(all_Hy,2));
end

%% Step 3: Cut region matrices to the same size
col_cut = 200;
sing_values = [];
for rr = 1:num_ranges
    num_cols = size(regions(rr).H_u, 2);
    if num_cols > col_cut       
        regions(rr).H_u   = regions(rr).H_u(:, 1:col_cut);
        regions(rr).H_y   = regions(rr).H_y(:, 1:col_cut);
    end
    fprintf("Region %d: Hankel-size: [%d  %d]\n", ...
    rr, size(regions(rr).H_u,1), size(regions(rr).H_u,2));
end
fprintf('\n=== HANKEL MATRIX CONSTRUCTION COMPLETE ===\n');

%% Step 4: Create overlapping region matrices
overlapping_regions = struct();
for rr = 1:1:num_ranges-1
    overlapping_regions(rr).H_u = [regions(rr).H_u, regions(rr+1).H_u];
    overlapping_regions(rr).H_y = [regions(rr).H_y, regions(rr+1).H_y];
    overlapping_regions(rr).low = regions(rr).range(1);
    overlapping_regions(rr).mid = regions(rr).range(2);
    overlapping_regions(rr).high = regions(rr+1).range(2);
    fprintf("Overlapping Region %d: Hankel-size: [%d  %d]\n", ...
    rr, size(overlapping_regions(rr).H_u,1), size(overlapping_regions(rr).H_u,2));
end

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
y_now = y0;

% Find the starting overlapping region:
for r_idx = 1:num_ranges-1
    lower = overlapping_regions(r_idx).low;
    upper = overlapping_regions(r_idx).high;
    if (r_idx < num_ranges &&  y_now >= lower && y_now < upper) || ...
       (r_idx == num_ranges && y_now >= lower && y_now <= upper)
        base_region = r_idx;
        break;
    end
end
%% Control loop:
for k = 1:ksim-T_fut
    % reference for the prediction horizon:
    y_reference = r(k:k+T_fut-1, :);
    %  Determine base region index for current output y_now 
    if (r_idx < num_ranges &&  y_now >= lower && y_now < upper) || ...
       (r_idx == num_ranges && y_now >= lower && y_now <= upper)
        % still in the same region: nothing to do
        base_region = r_idx; 
    else
        if y_now >= lower
            % left range trough the upper limit
            r_idx = r_idx + 1;
        else
            % left range trough the lower limit
            r_idx = r_idx - 1;
        end
        % out of range: find the new range
        lower = overlapping_regions(r_idx).low;
        upper = overlapping_regions(r_idx).high;
        if (r_idx < num_ranges &&  y_now >= lower && y_now < upper) || ...
           (r_idx == num_ranges && y_now >= lower && y_now <= upper)
            base_region = r_idx;
        end
    end
    if isnan(base_region)
        % if y_now somehow outside ranges, warn and fallback to region 1
        warning('k=%d: y_now=%.4f outside defined ranges. Falling back to region 1 Hankels.', k, y_now);
        H_u = overlapping_regions(1).H_u;
        H_y = overlapping_regions(1).H_y;
        H_idx(k) = 1;
    else
        if base_region < num_ranges
            % use determined base_region
            H_u = overlapping_regions(base_region).H_u;
            H_y = overlapping_regions(base_region).H_y;
            H_idx(k) = base_region; 
        else
            % use last overlapping region 
            H_u = overlapping_regions(num_ranges-1).H_u;
            H_y = overlapping_regions(num_ranges-1).H_y;
            H_idx(k) = num_ranges-1; 
        end
    end

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
    y_now = ydpc(k);
end
        
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

%% Save in RMSE vectors for plotting:
RMSE_SS(nn-1) = rmse_ss;
RMSE_trans(nn-1) = rmse_trans;
number_range(nn-1) = num_ranges;

%% Save results:
clear results;
%% NAME OF THIS EXPERIMENT 
experiment_name = ["GS_DeePC_3m_M" + num2str(num_ranges)];   % <-- CHANGE THIS FOR EACH RUN
results.name    = experiment_name;
results.ts      = ts;
results.tk      = tk(1:k);
results.r       = r(1:k);
results.ydpc    = ydpc(1:k);
results.H_idx   = H_idx(1:k);
results.rel_err = rel_err(1:k);
results.rel_err_full = rel_err_full(1:k);
results.rmse_ss = rmse_ss;
results.rmse_trans = rmse_trans;
results.mask = mask;
filename = "results_" + experiment_name + ".mat";
save(filename, "results");

fprintf('\nSaved results to: %s\n', filename);
fprintf('Steady-state RMSE: %.6f\n', rmse_ss);
fprintf('-------------------------------------------------\n');
end

function dx = unbaldisc(t,x,u) 
% dc motor with an unbalance disc from [2]
x1 = x(1);
x2 = x(2);

x1dot = x2;
x2dot = -1/0.5971*x2 - 0.07*9.8*0.42e-3/2.2e-4*sin(x1) + 15.3145/0.5971*u;

dx  = [x1dot; x2dot];
end


%% Plotting: 
% manually via Load_for_plot.m