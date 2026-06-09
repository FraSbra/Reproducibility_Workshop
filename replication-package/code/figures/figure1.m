%% BASELINE: IRF of EA interest rate to US MP shock
close all;

% -- Read data --
data = readtable("data_filtered.xlsx", "Sheet", "monthly");


eonia = data.EONIA; % ECB interest rate
ff4 = data.ff4; % Degasperi & Ricco (2021) MP shock (US)
hicp = log(data.HICP); % Log price level EA
ff_rate = data.FEDFUNDS; % Federal Funds Rate
IP_EA = log(data.IP_EA); % Industrial production EA
eu_dol = data.EU_DOL;
brent = log(data.Brent);

d_hicp = [NaN(1,1); diff(hicp)];
d_ff_rate = [NaN(1,1); diff(ff_rate)];
d_IP_EA = [NaN(1,1); diff(IP_EA)];
d_eonia = [NaN(1,1); diff(eonia)];
d_eu_dol = [NaN(1,1); diff(eu_dol)];
d_brent = [NaN(1,1); diff(brent)];

% -- Compute normalization constant --

y = d_ff_rate;
S = ff4;
H = 30; % Maximum horizon of IRF
hStart = 1; % Start LP at h = 0
lpMode = "level"; % Long-difference
seType = "ols";
c = 1; % Add constant
p = 12; % Maximum number of lags
X = lagmatrix([d_hicp d_brent], 0:p); % Set of controls
X = [X lagmatrix([d_ff_rate d_eonia d_IP_EA d_eu_dol], 1:p)];

res = lp_ols(y, X, S, H, hStart, lpMode, seType, c); % Run function
w = (1 / res.beta(1)) * 0.25;


%% IRF of euro interest rate
y = eonia;
S = d_ff_rate;
H = 30; % Maximum horizon of IRF
hStart = 0; % Start LP at h = 0
lpMode = "cum"; % Long-difference
c = 1; % Add constant
p = 12; % Maximum number of lags
Ziv = ff4; % Instrument 
X = lagmatrix([d_hicp d_brent], 0:p); % Set of controls
X = [X lagmatrix([d_ff_rate d_eonia d_IP_EA d_eu_dol], 1:p)];

res = lp_iv(y, X, S, Ziv, H, hStart, lpMode, c); % Run function

h     = res.h;
beta  = res.beta .* w;   % Normalized response
se    = res.se    .* w;  % Normalized standard errors

% Confidence bands
ub68 = beta + se;
lb68 = beta - se;
ub90 = beta + 1.645 * se;
lb90 = beta - 1.645 * se;

figure
hold on

x_fill = [h; flipud(h)];

% 90% band 
fill(x_fill, [ub90; flipud(lb90)], [0.60 0.75 1.00], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.4);

% 68% band (darker)
fill(x_fill, [ub68; flipud(lb68)], [0.20 0.50 0.90], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.4);

% IRF line and zero line
plot(h, beta,  'b-',  'LineWidth', 2);
yline(0,       'k--', 'LineWidth', 1, 'HandleVisibility', 'off');

xlabel('Horizon',                    'FontSize', 12)
ylabel('Response',                   'FontSize', 12)
title('IRF of EONIA to US MP shock', 'FontSize', 14)
legend('90% CI', '68% CI', 'Location', 'best')
grid on
box  on
hold off

% === SAVE ===
savePath = fullfile(saveResults_fig, 'eonia_irf.png');
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, savePath, '-dpng', '-r300');


