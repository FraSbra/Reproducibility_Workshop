% Horizon
H = 0:length(irf_fed)-1;

% Colors
redFed = [0.85 0.10 0.10];
blueECB = [0.10 0.30 0.85];

% Correct FED bands
lb68_fed = irf_fed - 1.00  * se_fed;
lb90_fed = irf_fed - 1.645 * se_fed;

figure;
hold on; box on; grid on;

%% ECB: blue filled confidence intervals

% 90% band
fill([H fliplr(H)], ...
     [ub90_ecb' fliplr(lb90_ecb')], ...
     blueECB, ...
     'FaceAlpha', 0.15, ...
     'EdgeColor', 'none');

% 68% band
fill([H fliplr(H)], ...
     [ub68_ecb' fliplr(lb68_ecb')], ...
     blueECB, ...
     'FaceAlpha', 0.30, ...
     'EdgeColor', 'none');

% ECB IRF
plot(H, irf_ecb, ...
     'Color', blueECB, ...
     'LineWidth', 2.2);

%% FED: red dashed intervals, not filled

% 90% dashed bounds
plot(H, ub90_fed, '--', 'Color', redFed, 'LineWidth', 1.4);
plot(H, lb90_fed, '--', 'Color', redFed, 'LineWidth', 1.4);

% 68% dashed bounds
plot(H, ub68_fed, ':', 'Color', redFed, 'LineWidth', 1.6);
plot(H, lb68_fed, ':', 'Color', redFed, 'LineWidth', 1.6);

% FED IRF
plot(H, irf_fed, ...
     'Color', redFed, ...
     'LineWidth', 2.2);

%% Zero line
yline(0, 'k-', 'LineWidth', 1);

%% Labels and legend
xlabel('Horizon');
ylabel('Response');
title('Downside Risk Response: FED vs ECB');

legend({'ECB 90% CI', 'ECB 68% CI', 'ECB IRF', ...
        'FED 90% CI', '', 'FED 68% CI', '', 'FED IRF'}, ...
        'Location', 'best');

hold off;

% === SAVE ===
savePath = fullfile(saveResults_fig, 'comparison_downside_risk_FED_vs_ECB.png');
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, savePath, '-dpng', '-r300');
