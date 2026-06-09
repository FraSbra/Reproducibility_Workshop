%% Figure A2: IRF of downside risk to INFO shock
clear;
close all;
clc;

addpath(genpath("../../Excel")); % Folder containing excel files
addpath(genpath("../routines")); % Folder contains all functions
appendixDir = fileparts(mfilename("fullpath"));
saveResults = fullfile(appendixDir, "results_appendix");


%% Import
other = readtable("data_filtered.xlsx", "Sheet", "other");
gdp_panel = readtable("data_filtered.xlsx", "Sheet", "GDP_panel");
load("pred_q.mat");
load("yy.mat");

CountryNames = gdp_panel.Properties.VariableNames(2:end);
fedfunds = other.FEDFUNDS;
eonia = other.EONIA;
eu_dol = other.DEXUSEU;
ciss = other.CISS;
us_gdp = log(other.GDPC1);
info = other.INFO;

%% Setup LP
settingsLP.maxH = 20; % Maximum number of horizons
settingsLP.MaxLPLagsOwn = 4; % Maximum number of lags of dependent variable
settingsLP.MaxLPLagsOther = 4; % Maximum number of lags of other controls
interaction = 'No_Interactions';

%% Common inputs

No_Countries = size(CountryNames,2); % Number of countries
Tt = size(gdp_panel,1);

Data.id = repmat(1:No_Countries,Tt,1); clear Tt
Data.quarter = repmat((1:size(gdp_panel,1))',1,No_Countries);

%% Structure

% Dependent
yy_d = squeeze(pred_q(:,10,:) - pred_q(:,1,:));

% Shock and instrumented variable
Instrument = info;
XX = fedfunds;

% Controls
Controls = cat(3, ...
    MatrixDiff(repmat(XX,1,No_Countries)), ...
    MatrixDiff(repmat(100 .* log(eu_dol),1,No_Countries)), ...
    MatrixDiff(repmat(eonia,1,No_Countries)), ...
    MatrixDiff(repmat(ciss,1,No_Countries)), ...
    MatrixDiff(repmat(us_gdp,1,No_Countries)) ...
);

%% LOCAL PROJECTIONS
Collect_LP = Panel_LP_IV(yy_d, XX, Controls, [], Instrument, Data.id, Data.quarter, ...
    settingsLP, interaction);

Collect.Coeff(:) = Collect_LP.Coeff .* 0.25 / Collect_LP.delta_hat; % IRFs
Collect.StdErr(:) = Collect_LP.StdErr .* 0.25 / Collect_LP.delta_hat; % Standard errors

%% Plot IRFs
h = (0:settingsLP.maxH)';
irf = Collect.Coeff(:);
se  = Collect.StdErr(:);

% 68% bands (~1 std err)
ub68 = irf + 1.00 * se;
lb68 = irf - 1.00 * se;

% 90% bands (~1.645 std err)
ub90 = irf + 1.645 * se;
lb90 = irf - 1.645 * se;

figure;
hold on;

x_fill = [h; flipud(h)];

% 90% band (lighter blue, plotted first so 68% sits on top)
fill(x_fill, [ub90; flipud(lb90)], [0.60 0.75 1.00], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.35);

% 68% band (darker blue)
fill(x_fill, [ub68; flipud(lb68)], [0.20 0.50 0.90], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.40);

plot(h, irf, 'b-', 'LineWidth', 2, 'HandleVisibility', 'off');
yline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');

xlabel('Horizon');
ylabel('Response');
title('IRF of downside risk to INFO shock');
legend('90% CI', '68% CI', 'Location', 'best');
grid on;
hold off;

% === SAVE FIGURE AS PNG ===
savePath = fullfile(saveResults, 'linear_downside_risk_irf_INFO.png');
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, savePath, '-dpng', '-r300');


%% Comparison with average effect

Collect_LP_AVG = Panel_LP_IV(yy, XX, Controls, [], Instrument, Data.id, Data.quarter, ...
    settingsLP, interaction);

Collect_AVG.Coeff(:) =  - Collect_LP_AVG.Coeff .* 0.25 / Collect_LP_AVG.delta_hat; % IRFs --> minus in front for comparison
Collect_AVG.StdErr(:) = Collect_LP_AVG.StdErr .* 0.25 / Collect_LP_AVG.delta_hat; % Standard errors

%% Comparison plot: Downside Risk vs Average Effect
h = (0:settingsLP.maxH)';

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
%title('IRF to INFO Shock');
legend('Location', 'best');
grid on;
hold off;

% === SAVE ===
savePath = fullfile(saveResults, 'comparison_downside_risk_vs_avg_irf_INFO.png');
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, savePath, '-dpng', '-r300');
