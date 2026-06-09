%% Figure A1: IRF of US downside risk to US MP shock
close all;
clc;

addpath(genpath("../../excel"));
addpath(genpath("../routines"));

appendixDir = fileparts(mfilename("fullpath"));
saveResults = fullfile(appendixDir, "results_appendix");

%% Import

other = readtable("data_filtered.xlsx", "Sheet", "other");
load(fullfile(saveResults, "pred_q_usa.mat"));

% Create mask starting from January 1999
mask = other.observation_date >= datetime(1999,1,1);

date = other.observation_date(mask);
pred_q = pred_q(mask,:,:);

fedfunds = other.FEDFUNDS(mask);
eonia = other.EONIA(mask);
eu_dol = other.DEXUSEU(mask);
eu_shock = other.eu_shock(mask);
ff4 = other.ff4(mask);
gdpusa = 100*log(other.GDPC1(mask));

%% Setup LP
H = 20; % Maximum horizon of IRF
hStart = 0; % Start LP at h = 0
lpMode = "cum"; % Long-difference
c = 1; % Add constant
p = 4; % Maximum number of lags

%% Structure

% Dependent
yy_d = pred_q(:,10) - pred_q(:,1);

% Shock and instrumented variable
Instrument = ff4;
XX = MatrixDiff(fedfunds);

% Controls
Controls = lagmatrix([ ...
    MatrixDiff(fedfunds), ...
    ff4 ...
    MatrixDiff(gdpusa) ...
], 1:p);

%% LOCAL PROJECTIONS
res = lp_iv(yy_d, Controls, XX, Instrument, H, hStart, lpMode, c);

% First stage for normalization
y = MatrixDiff(fedfunds);
S = ff4;
res_fs = lp_ols(y, Controls, S, H, 1, "level", "ols", c);
w = 0.25 / (res_fs.beta(1));

Collect.Coeff(:) = res.beta .* w; % IRFs
Collect.StdErr(:) = res.se .* w; % Standard errors

%% Plot IRFs
h = res.h;
irf_fed = Collect.Coeff(:);
se_fed  = Collect.StdErr(:);

% 68% bands (~1 std err)
ub68_fed = irf_fed + 1.00 * se_fed;
lb68_fed = irf_fed - 1.00 * se_fed;

% 90% bands (~1.645 std err)
ub90_fed = irf_fed + 1.645 * se_fed;
lb90_fed = irf_fed - 1.645 * se_fed;

figure;
hold on;

x_fill = [h; flipud(h)];

% 90% band (lighter blue, plotted first so 68% sits on top)
fill(x_fill, [ub90_fed; flipud(lb90_fed)], [0.60 0.75 1.00], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.35);

% 68% band (darker blue)
fill(x_fill, [ub68_fed; flipud(lb68_fed)], [0.20 0.50 0.90], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.40);

plot(h, irf_fed, 'b-', 'LineWidth', 2, 'HandleVisibility', 'off');
yline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');

xlabel('Horizon');
ylabel('Response');
title('IRF of US downside risk to EU MP shock');
legend('90% CI', '68% CI', 'Location', 'best');
grid on;
hold off;

% === SAVE FIGURE AS PNG ===
savePath = fullfile(saveResults, 'linear_downside_risk_irf_USA.png');
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, savePath, '-dpng', '-r300');


%% Comparison with average effect


yy = yy(mask);
res_avg = lp_iv(gdpusa, Controls, XX, Instrument, H, hStart, lpMode, c);

Collect_AVG.Coeff(:) = - res_avg.beta .* w; % IRFs --> minus in front for comparison
Collect_AVG.StdErr(:) = res_avg.se .* w; % Standard errors

%% Comparison plot: Downside Risk vs Average Effect
h = res.h;

% --- Downside risk ---
irf1 = Collect.Coeff(:);
se1  = Collect.StdErr(:);
ub68_1 = irf1 + 1.00 * se1;  lb68_1 = irf1 - 1.00 * se1;
ub90_1 = irf1 + 1.645 * se1; lb90_1 = irf1 - 1.645 * se1;

% --- Average effect ---
irf2 = Collect_AVG.Coeff(:);
se2  = Collect_AVG.StdErr(:);
ub68_2 = irf2 + 1.00 * se2;  lb68_2 = irf2 - 1.00 * se2;
ub90_2 = irf2 + 1.645 * se2; lb90_2 = irf2 - 1.645 * se2;

x_fill = [h; flipud(h)];

figure;
hold on;

% 90% bands
fill(x_fill, [ub90_1; flipud(lb90_1)], [0.60 0.75 1.00], 'EdgeColor', 'none', 'FaceAlpha', 0.35, 'HandleVisibility', 'off');
fill(x_fill, [ub90_2; flipud(lb90_2)], [1.00 0.75 0.60], 'EdgeColor', 'none', 'FaceAlpha', 0.35, 'HandleVisibility', 'off');

% 68% bands
fill(x_fill, [ub68_1; flipud(lb68_1)], [0.20 0.50 0.90], 'EdgeColor', 'none', 'FaceAlpha', 0.40, 'HandleVisibility', 'off');
fill(x_fill, [ub68_2; flipud(lb68_2)], [0.90 0.40 0.20], 'EdgeColor', 'none', 'FaceAlpha', 0.40, 'HandleVisibility', 'off');

% IRF lines
plot(h, irf1, 'b-', 'LineWidth', 2, 'DisplayName', 'Downside Risk');
plot(h, irf2, 'r-', 'LineWidth', 2, 'DisplayName', 'Average Effect');

yline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');

xlabel('Horizon');
ylabel('Response');
title('IRF to US MP Shock');
legend('Location', 'best');
grid on;
hold off;

% === SAVE ===
savePath = fullfile(saveResults, 'comparison_downside_risk_vs_avg_irf_USA.png');
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, savePath, '-dpng', '-r300');
