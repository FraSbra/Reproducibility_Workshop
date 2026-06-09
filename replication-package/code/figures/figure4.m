%% Figure 4: Counterfactual --> no financial response
close all

FWDControls = repmat(ciss,1,No_Countries);


% --- Counterfactual ---
Collect_LPCounter=Panel_LP_IV_withFWDControls(yy_d,XX,FWDControls,Controls,...
    [],Instrument,Data.id,Data.quarter,settingsLP, interaction);
    
Counter.Coeff(:) = Collect_LPCounter.Coeff .* 0.25 / Collect_LPCounter.delta_hat; % IRFs
Counter.StdErr(:) = Collect_LPCounter.StdErr .* 0.25 / Collect_LPCounter.delta_hat; % Standard errors

%% Plot IRFs
h = (0:settingsLP.maxH)';
irf = Collect.Coeff(:);
se  = Collect.StdErr(:);

% Counterfactual
irf_c = Counter.Coeff(:);
se_c  = Counter.StdErr(:);

% 68% bands (~1 std err)
ub68 = irf + 1.00 * se;
lb68 = irf - 1.00 * se;

% 90% bands (~1.645 std err)
ub90 = irf + 1.645 * se;
lb90 = irf - 1.645 * se;

% Counterfactual bands
ub68_c = irf_c + 1.00 * se_c;
lb68_c = irf_c - 1.00 * se_c;
ub90_c = irf_c + 1.645 * se_c;
lb90_c = irf_c - 1.645 * se_c;

figure;
hold on;

x_fill = [h; flipud(h)];

% 90% band baseline (lighter blue)
fill(x_fill, [ub90; flipud(lb90)], [0.60 0.75 1.00], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.35, 'HandleVisibility', 'off');

% 68% band baseline (darker blue)
fill(x_fill, [ub68; flipud(lb68)], [0.20 0.50 0.90], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.40, 'HandleVisibility', 'off');

% Baseline IRF
plot(h, irf, 'b-', 'LineWidth', 2);

% Counterfactual IRF
plot(h, irf_c, 'k-', 'LineWidth', 2);

% Counterfactual 68% bands
plot(h, ub68_c, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
plot(h, lb68_c, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');

% Counterfactual 90% bands
plot(h, ub90_c, 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');
plot(h, lb90_c, 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');

yline(0, 'k-', 'LineWidth', 0.5, 'HandleVisibility', 'off');

xlabel('Horizon');
ylabel('Response');
title('IRF of downside risk to US MP shock');
legend('Baseline', 'Counterfactual', 'Location', 'best');
grid on;
hold off;

% === SAVE FIGURE AS PNG ===
savePath = fullfile(saveResults_fig, 'counter_no_CISS.png');
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, savePath, '-dpng', '-r300');
