dbstop if error
clear; clc;

% load retweeting probability
load('./example.mat');
x = D1(:, 1);
y = D1(:, 2);
gamma = 0.24;


%% Fitting Method1: only estimates omega
alpha = y(1) / (1-gamma);
[Goodness, Paras, ci_1, ci_2] = fit_BetaX_Method1(alpha, gamma, x, y);
RMSE = Goodness(1); 
R2ADJ = Goodness(3);

figure; hold on
set(gcf, 'Position', [100 200 550 450]);
set(gcf, 'Color', 'White');

omega = Paras;
y_fit = alpha * x .* ((1-gamma) .^ (x.^omega));
fill([x', fliplr(x')], [ci_2(:,1)', fliplr(ci_2(:,2)')], [.9 .9 .9], 'EdgeColor', [.9 .9 .9]);
h1 = plot(x, y, 'bo', 'LineWidth', 1.5, 'MarkerSize', 8); 
h2 = plot(x, y_fit, 'r-', 'LineWidth', 1.5);
text(2, 3e-3, {strcat('\alpha=',num2str(round(alpha,4)))}, 'FontName', 'Helvetica', 'FontSize', 16);
text(2, 2e-3, {strcat('\omega=',num2str(round(omega,2)))}, 'FontName', 'Helvetica', 'FontSize', 16);
xlabel('Number of Exposures, x', 'FontName','Helvetica');
ylabel('Retweeting Probability, \beta(x)', 'FontName','Helvetica');
title('Fitting Method1', 'FontName','Helvetica');
legend([h1 h2], {'Empirical', 'Fitting'}, 'FontName','Helvetica', 'FontSize', 16, 'Location', 'NE', 'Box', 'Off');
set(gca, 'FontName','Helvetica', 'FontSize', 16, 'Box', 'On', 'LineWidth',1);
set(gca, 'XTickLabelMode', 'auto', 'YTickLabelMode', 'auto');
ax = gca; 
ax.YAxis.Exponent = -3;
xlim([0 7]); ylim([0 1e-2]);


%% Fitting Method2: estimates both alpha and omega
[Goodness, Paras, ci_1, ci_2] = fit_BetaX_Method2(gamma, x, y);
RMSE = Goodness(1); 
R2ADJ = Goodness(3);

figure; hold on
set(gcf, 'Position', [300 200 550 450]);
set(gcf, 'Color', 'White');

alpha = Paras(1); omega = Paras(2);
y_fit = alpha * x .* ((1-gamma) .^ (x.^omega));
fill([x', fliplr(x')], [ci_2(:,1)', fliplr(ci_2(:,2)')], [.9 .9 .9], 'EdgeColor', [.9 .9 .9]);
h1 = plot(x, y, 'bo', 'LineWidth', 1.5, 'MarkerSize', 8); 
h2 = plot(x,y_fit,'r-', 'LineWidth', 1.5);
text(2, 3e-3, {strcat('\alpha=',num2str(round(alpha,4)))}, 'FontName', 'Helvetica', 'FontSize', 16);
text(2, 2e-3, {strcat('\omega=',num2str(round(omega,2)))}, 'FontName', 'Helvetica', 'FontSize', 16);
xlabel('Number of Exposures, x', 'FontName','Helvetica');
ylabel('Retweeting Probability, \beta(x)', 'FontName','Helvetica');
title('Fitting Method2', 'FontName','Helvetica');
legend([h1 h2], {'Empirical', 'Fitting'}, 'FontName','Helvetica', 'FontSize', 16, 'Location', 'NE', 'Box', 'Off');
set(gca, 'FontName','Helvetica', 'FontSize', 16, 'Box', 'On', 'LineWidth',1);
set(gca, 'XTickLabelMode', 'auto', 'YTickLabelMode', 'auto');
ax = gca; 
ax.YAxis.Exponent = -3;
xlim([0 7]); ylim([0 1e-2]);


%%  Fitting Method3: estimates alpha, gamma, and omega all together
BetaX = load('./example.mat');
gamma = 0.1 : 0.02 : 0.8;
s = length(fieldnames(BetaX));
RMSE = zeros(s, length(gamma)); 
R2ADJ = zeros(s, length(gamma));
alpha_est = zeros(s, length(gamma));
omega_est = zeros(s, length(gamma));

for i = 1 : s
    disp(strcat('message-', num2str(i)));
    D = getfield(BetaX, strcat('D', num2str(i)));
    x = D(:, 1);
    y = D(:, 2);
    [R1, ~, R2, ~] = grid_search_optimal_gamma(gamma, x, y);
    RMSE(i, :) = R1(:, 1)';
    R2ADJ(i, :) = R1(:, 2)';
    alpha_est(i, :) = R2(:, 1)';
    omega_est(i, :) = R2(:, 2)';
end

RMSE_avg = mean(RMSE, 1);
R2ADJ_avg = mean(R2ADJ, 1);
idx = find(RMSE_avg==min(RMSE_avg));

figure; hold on
set(gcf, 'Position', [500 200 550 450]);
set(gcf, 'Color', 'White');
plot(gamma, RMSE_avg, 'b-o', 'LineWidth', 1);
hold on
plot([gamma(idx) gamma(idx)], [0 1], 'k--', 'LineWidth', 1);
xlabel('Proportion of Common Neighbors, \gamma', 'FontName','Helvetica');
ylabel('RMSE', 'FontName','Helvetica');
title('Fitting Method3', 'FontName','Helvetica');
set(gca, 'FontName','Helvetica', 'FontSize', 16, 'Box', 'On', 'LineWidth',1);
set(gca, 'XTickLabelMode', 'auto', 'YTickLabelMode', 'auto');
ax = gca; 
ax.YAxis.Exponent = -4;
xlim([0.08 0.82]); ylim([3e-4 12e-4]);


%% grid search optimal gamma
function [GoF_AVG, GoF_STD, Paras_AVG, Paras_STD] = grid_search_optimal_gamma(gamma, x, y)
    maxIter = 10; % maximum number of fitting
    
    GoF_AVG = zeros(length(gamma), 2);
    GoF_STD = zeros(length(gamma), 2);
    Paras_AVG = zeros(length(gamma), 2);
    Paras_STD = zeros(length(gamma), 2);
    
    disp('trying gamma:');
    for i = 1 : length(gamma)
        disp(gamma(i));
        temp_GoF = zeros(maxIter, 2); 
        temp_Paras = zeros(maxIter, 2);
        for j = 1 : maxIter
            [gof, Paras, ~, ~] = fit_BetaX_Method2(gamma(i), x, y);
            temp_GoF(j, :) = [gof(1), gof(3)];
            temp_Paras(j, :) = [Paras(1), Paras(2)];
        end
        GoF_AVG(i, :) = mean(temp_GoF);
        GoF_STD(i, :) = std(temp_GoF, 1);
        Paras_AVG(i, :) = mean(temp_Paras);
        Paras_STD(i, :) = std(temp_Paras, 1);
    end
end