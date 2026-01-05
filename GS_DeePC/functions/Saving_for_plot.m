%% File Name: Saving_for_plot.m
% Author: Sebastian Zieglmeier
% Date last updated: 30.11.2025
% Description: Manually save results in mat file
%
% Sources:
% [1] - Sebastian Zieglmeier, et al., "Gain-Scheduling Data-Enabled 
%       Predictive Control for Nonlinear Systems with Linearized Operating 
%       Regions", 
%       https://doi.org/10.48550/arXiv.2512.02797
%  
%  Save function for manually saving the results to reload and plot them
%  later

clear results;

%% === NAME OF THIS EXPERIMENT ===================================
experiment_name = "GS_DeePC_no_overlap_no_dwell_time_M8";  
% =================================================================

%% === Compute RMSE only inside steady-state intervals =============

% mask: where relative error is non-zero (steady-state intervals)
mask = rel_err ~= 0;

%% === Build results struct ========================================

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

%% === Save file ====================================================

filename = "results_" + experiment_name + ".mat";
save(filename, "results");

fprintf('\nSaved results to: %s\n', filename);
fprintf('Steady-state RMSE: %.6f\n', rmse_ss);
fprintf('-------------------------------------------------\n');
