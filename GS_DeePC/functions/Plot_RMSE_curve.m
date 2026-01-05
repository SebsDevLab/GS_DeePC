%% File Name: Plot_RMSE_curve.m
% Author: Sebastian Zieglmeier
% Date last updated: 30.11.2025
% Description: Plot RMSE curve 
%
% Sources:
% [1] - Sebastian Zieglmeier, et al., "Gain-Scheduling Data-Enabled 
%       Predictive Control for Nonlinear Systems with Linearized Operating 
%       Regions", 
%       https://doi.org/10.48550/arXiv.2512.02797
%  
%  Load and plot the RMSEs of the looped GS_DeePC

RMSE_SS = [];
RMSE_T = [];
for nn=2:1:50
    
    result_string = ['results_GS_DeePC_M'  num2str(nn)  '.mat'];
    try
        load(result_string);
        RMSE_SS(nn-1) = results.rmse_ss;
        RMSE_T(nn-1) = results.rmse_trans;
        
    catch
        warning(["Could not load: " + result_string]);
    end
end

% --- Publication-quality plot of RMSE metrics ---

figure('Color','w');                % white background
hold on;

% --- Left axis: RMSE_SS ---
yyaxis left
p1 = plot(RMSE_SS, 'LineWidth', 1.8);
p1.Color = [0 0.447 0.741];         % MATLAB default blue (good visibility)

xlabel('$C_n$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$RMSE_{SS}$ [rad]', 'Interpreter','latex', 'FontSize', 12)

% --- Right axis: RMSE_T ---
yyaxis right
p2 = plot(RMSE_T, 'LineWidth', 1.8);
p2.Color = [0.850 0.325 0.098];     % MATLAB default orange

ylabel('$RMSE_{T}$ [rad]', 'Interpreter','latex', 'FontSize', 12)

% --- Style tweaks ---
grid on
box on
set(gca, 'FontSize', 11)

% --- High-quality export ---
print(gcf, 'RMSE_plot.png', '-dpng', '-r600');   % 600 dpi PNG

