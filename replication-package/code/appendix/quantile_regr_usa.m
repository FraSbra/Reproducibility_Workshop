%% Compute predicted quantiles for US GDP

%% Notation:
%   V1: NFCI MIDAS & GDP lag
%   _AR: augmented with autoregressive term

clear;
close all;
clc;

addpath(genpath("../../excel"));
addpath(genpath("../routines"));
appendixDir = fileparts(mfilename("fullpath"));
saveResults = fullfile(appendixDir, "results_appendix");



other = readtable("data_filtered.xlsx", "Sheet", "other");
monthly = readtable("data_filtered.xlsx", "Sheet", "monthly");

date = other.observation_date;
monthly_date = monthly.observation_date;
nfci_m = monthly.NFCI;
gdp = 100* MatrixDiff(log(other.GDPC1));

gdp_lag = lagmatrix(gdp,1);
nfci_midas = NaN(size(gdp,1), 3);

for m = 1:3 % Loop over the three monthly lags within the previous quarter
    [~, idx] = ismember(date - calmonths(4-m), monthly_date); % Match each quarterly date with the corresponding previous monthly date
    nfci_midas(idx > 0,m) = nfci_m(idx(idx > 0)); % Store NFCI monthly values where the monthly date exists
end

nfci_midas_lag = lagmatrix(nfci_midas,1); %trust me -> perdiamo un osservazione da growth rate gdp (per quello non andava prima)

quantiles = 0.05:0.05:0.95;

yy = NaN(size(gdp));

for t=1:size(gdp,1) - 3
    yy(t) = mean(gdp(t:t+3));
end

RHS_V1 = [ones(size(gdp,1),1) nfci_midas_lag gdp_lag];

pred_q_V1 = NaN(size(gdp,1), length(quantiles));
pred_q_V1_AR = NaN(size(gdp,1), length(quantiles));

% --- V1 ---
i = 1;
beta_hat = NaN(size(RHS_V1,2), length(quantiles));
for q = quantiles
    RHS = RHS_V1;
    beta_hat(:,i) = smooth_quantile_reg(yy, RHS, q);
    [~, ~, pred_q_V1_AR(:,i)] = smooth_quantile_reg_ar(yy, RHS, q);
    i = i + 1;
end
pred_q_V1 = RHS * beta_hat;

pred_q = pred_q_V1;
pred_q_AR = pred_q_V1_AR;
downside = -(pred_q_AR(:,10) - pred_q_AR(:,2));
upside   =   pred_q_AR(:,18) - pred_q_AR(:,10);

save(fullfile(saveResults, "pred_q_usa.mat"), "pred_q", "pred_q_AR", ...
    "pred_q_V1", "pred_q_V1_AR", "yy", "downside", "upside", ...
    "date", "gdp", "nfci_midas");

%% Figure: Downside and upside plots

figure;

colors = [0.2 0.4 0.8;   % blue for downside
          0.8 0.2 0.2];  % red for upside

recessions = [datetime(1974,10,1), datetime(1975,3,31);
              datetime(1980,1,1),  datetime(1980,9,30);
              datetime(1982,1,1),  datetime(1982,12,31);
              datetime(1992,1,1),  datetime(1993,3,31);
              datetime(2008,1,1),  datetime(2009,6,30);
              datetime(2011,7,1),  datetime(2013,3,31);
              datetime(2020,1,1),  datetime(2020,6,30)];

ymin = min([downside; upside]);
ymax = max([downside; upside]);

% Shade recession bands first
for r = 1:size(recessions, 1)
    fill([recessions(r,1), recessions(r,2), recessions(r,2), recessions(r,1)], ...
         [ymin, ymin, ymax, ymax], ...
         [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    hold on;
end

h1 = plot(date, downside, 'Color', colors(1,:), 'LineWidth', 1.2);
h2 = plot(date, upside,   'Color', colors(2,:), 'LineWidth', 1.2);
hold off;

title('United States', 'FontSize', 9);
xlim([datetime(1980,1,1), date(end)]);
ylim([ymin, ymax]);
grid on;
box off;

legend([h1, h2], {'Downside uncertainty', 'Upside uncertainty'}, ...
    'Orientation', 'horizontal', 'FontSize', 9, 'Location', 'southoutside');

% === SAVE FIGURE AS PNG ===
savePath = fullfile(saveResults, 'down_up_USA.png');
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, savePath, '-dpng', '-r300');
