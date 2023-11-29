dbstop if error
clear; clc;

gamma = 0.24;

%% plot the first ten subfigures of FIGS34
for mid = 1 : 10
    
    load(strcat('./FIGS34_subfigure_', num2str(mid), '.mat'));
    x = D(:, 1);
    y = D(:, 2);
    
    [Goodness, Paras, ci_1, ci_2] = fit_BetaX_Method2(gamma, x, y);
    RMSE = Goodness(1); R2ADJ = Goodness(3);
    
    figure; hold on
    set(gcf,'Color','White');

    alpha = Paras(1); omega = Paras(2);
    y_fit = alpha * x .* ((1-gamma).^(x.^omega));
    fill([x', fliplr(x')], [ci_2(:,1)', fliplr(ci_2(:,2)')], [.9 .9 .9], 'EdgeColor', [.9 .9 .9]);
    h1 = plot(x, y, 'bo', 'LineWidth', 1.5, 'MarkerSize', 8); 
    h2 = plot(x, y_fit, 'r-', 'LineWidth', 1.5);
    text(length(x)-2, 1.2*max(y), {strcat('\alpha=',num2str(round(alpha,4)))},'FontName','Helvetica','FontSize',16);
    text(length(x)-2, max(y), {strcat('\omega=',num2str(round(omega,2)))},'FontName','Helvetica','FontSize',16);
    xlabel('Number of Exposures, x', 'FontName','Helvetica');
    ylabel('Retweeting Probability, \beta(x)', 'FontName','Helvetica');
    title(strcat('subfigure-', num2str(mid)), 'FontName','Helvetica');
    set(gca, 'FontName', 'Helvetica', 'FontSize', 16, 'Box', 'On', 'LineWidth', 1);
    set(gca,'XTickLabelMode', 'auto', 'YTickLabelMode', 'auto');
    ax = gca; ax.YAxis.Exponent = -3;
    xlim([0 length(x)+1]); ylim([0 1.5*max(y)]);

end