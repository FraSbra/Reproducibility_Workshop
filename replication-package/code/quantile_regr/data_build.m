%% Compute predicted quantiles 

%% Notation:
%   V1: NFCI & GDP
%   V2: NFCI, CISS & GDP
%   V3: CISS & GDP
%   V4: CISS MIDAS, NFCI MIDAS & GDP
%   V5: NFCI, CISS, debt/GDP & GDP
%   V6: NFCI MIDAS, CISS & GDP
%   V7: CISS MIDAS & GDP
%   _AR: augmented with autoregressive term

gdp_panel = readtable("data_filtered.xlsx", "Sheet", "GDP_panel");
other = readtable("data_filtered", "Sheet", "other");
monthly = readtable("data_filtered", "Sheet", "monthly");
debt_panel = readtable("data_filtered", "Sheet", "gov");
clc

date = other.observation_date;
monthly_date = monthly.observation_date;
countryNames = gdp_panel.Properties.VariableNames(2:end);
nfci = other.NFCI;
ciss = other.CISS;
nfci_m = monthly.NFCI;
ciss_m = monthly.CISS;

gdp = table2array(gdp_panel(:,2:end));
gdp_lag = lagmatrix(gdp,1);
debt = table2array(debt_panel(:,2:end));
debt_lag = lagmatrix([NaN(1,size(debt,2)); 100*diff(log(debt))],1);


nfci_lag = lagmatrix(nfci,1);
ciss_lag = lagmatrix(ciss,1);
nfci_midas = NaN(size(gdp,1), 3);
ciss_midas = NaN(size(gdp,1), 3);

for m = 1:3 % Loop over the three monthly lags within the previous quarter
    [~, idx] = ismember(date - calmonths(4-m), monthly_date); % Match each quarterly date with the corresponding previous monthly date
    nfci_midas(idx > 0,m) = nfci_m(idx(idx > 0)); % Store NFCI monthly values where the monthly date exists
    ciss_midas(idx > 0,m) = ciss_m(idx(idx > 0)); % Store CISS monthly values where the monthly date exists
end

quantiles = 0.05:0.05:0.95;

yy = NaN(size(gdp));

for i=1:size(gdp,2)
    for t=1:size(gdp,1) - 3
        yy(t,i) = sum(gdp(t:t+3,i));
    end
end

% Save variable (sum of gdp within year) to compute average effects of shocks
save("../../excel/processed/yy.mat", "yy"); % time x country

RHS_V1 = [ones(size(gdp,1),1) nfci_lag];
RHS_V2 = [ones(size(gdp,1),1) nfci_lag ciss_lag];
RHS_V3 = [ones(size(gdp,1),1) ciss_lag];
RHS_V4 = [ones(size(gdp,1),1) ciss_midas nfci_midas];
RHS_V6 = [ones(size(gdp,1),1) nfci_midas ciss_lag];
RHS_V7 = [ones(size(gdp,1),1) ciss_midas];

pred_q_V1 = NaN(size(gdp,1), length(quantiles), size(gdp,2));
pred_q_V2 = NaN(size(gdp,1), length(quantiles), size(gdp,2));
pred_q_V3 = NaN(size(gdp,1), length(quantiles), size(gdp,2));
pred_q_V4 = NaN(size(gdp,1), length(quantiles), size(gdp,2));
pred_q_V5 = NaN(size(gdp,1), length(quantiles), size(gdp,2));
pred_q_V6 = NaN(size(gdp,1), length(quantiles), size(gdp,2));
pred_q_V7 = NaN(size(gdp,1), length(quantiles), size(gdp,2));

pred_q_V1_AR = NaN(size(gdp,1), length(quantiles), size(gdp,2));
pred_q_V2_AR = NaN(size(gdp,1), length(quantiles), size(gdp,2));
pred_q_V3_AR = NaN(size(gdp,1), length(quantiles), size(gdp,2));
pred_q_V4_AR = NaN(size(gdp,1), length(quantiles), size(gdp,2));
pred_q_V5_AR = NaN(size(gdp,1), length(quantiles), size(gdp,2));
pred_q_V6_AR = NaN(size(gdp,1), length(quantiles), size(gdp,2));
pred_q_V7_AR = NaN(size(gdp,1), length(quantiles), size(gdp,2));

% --- V1 ---
for j =1:size(gdp,2)
    i = 1;
    beta_hat = NaN(size(RHS_V1,2) + 1, length(quantiles));
    RHS = [RHS_V1 gdp_lag(:,j)];
    for q = quantiles
        beta_hat(:,i) = smooth_quantile_reg(yy(:,j), RHS, q);
        [~, ~, pred_q_V1_AR(:,i,j)] = smooth_quantile_reg_ar(yy(:,j), RHS, q);
        i = i + 1;
    end
    pred_q_V1(:,:,j) = RHS * beta_hat;
end

% --- V2 ---
for j =1:size(gdp,2)
    i = 1;
    beta_hat = NaN(size(RHS_V2,2) + 1, length(quantiles));
    for q = quantiles
        RHS = [RHS_V2 gdp_lag(:,j)];
        beta_hat(:,i) = smooth_quantile_reg(yy(:,j), RHS, q);
        [~, ~, pred_q_V2_AR(:,i,j)] = smooth_quantile_reg_ar(yy(:,j), RHS, q);
        i = i + 1;
    end
    pred_q_V2(:,:,j) = RHS * beta_hat;
end

% --- V3 ---
for j =1:size(gdp,2)
    i = 1;
    beta_hat = NaN(size(RHS_V3,2) + 1, length(quantiles));
    RHS = [RHS_V3 gdp_lag(:,j)];
    for q = quantiles
        beta_hat(:,i) = smooth_quantile_reg(yy(:,j), RHS, q);
        [~, ~, pred_q_V3_AR(:,i,j)] = smooth_quantile_reg_ar(yy(:,j), RHS, q);
        i = i + 1;
    end
    pred_q_V3(:,:,j) = RHS * beta_hat;
end

% --- V4 ---
for j =1:size(gdp,2)
    i = 1;
    beta_hat = NaN(size(RHS_V4,2) + 1, length(quantiles));
    RHS = [RHS_V4 gdp_lag(:,j)];
    for q = quantiles
        beta_hat(:,i) = smooth_quantile_reg(yy(:,j), RHS, q);
        [~, ~, pred_q_V4_AR(:,i,j)] = smooth_quantile_reg_ar(yy(:,j), RHS, q);
        i = i + 1;
    end
    pred_q_V4(:,:,j) = RHS * beta_hat;
end

% --- V5 ---
for j =1:size(gdp,2)
    i = 1;
    beta_hat = NaN(size(RHS_V2,2) + 2, length(quantiles));
    RHS = [RHS_V2 debt_lag(:,j) gdp_lag(:,j)];
    for q = quantiles  
        beta_hat(:,i) = smooth_quantile_reg(yy(:,j), RHS, q);
        [~, ~, pred_q_V5_AR(:,i,j)] = smooth_quantile_reg_ar(yy(:,j), RHS, q);
        i = i + 1;
    end
    pred_q_V5(:,:,j) = RHS * beta_hat;
end

% --- V6 ---
for j =1:size(gdp,2)
    i = 1;
    beta_hat = NaN(size(RHS_V6,2) + 1, length(quantiles));
    RHS = [RHS_V6 gdp_lag(:,j)];
    for q = quantiles
        beta_hat(:,i) = smooth_quantile_reg(yy(:,j), RHS, q);
        [~, ~, pred_q_V6_AR(:,i,j)] = smooth_quantile_reg_ar(yy(:,j), RHS, q);
        i = i + 1;
    end
    pred_q_V6(:,:,j) = RHS * beta_hat;
end

% --- V7 ---
for j =1:size(gdp,2)
    i = 1;
    beta_hat = NaN(size(RHS_V7,2) + 1, length(quantiles));
    RHS = [RHS_V7 gdp_lag(:,j)];
    for q = quantiles
        beta_hat(:,i) = smooth_quantile_reg(yy(:,j), RHS, q);
        [~, ~, pred_q_V7_AR(:,i,j)] = smooth_quantile_reg_ar(yy(:,j), RHS, q);
        i = i + 1;
    end
    pred_q_V7(:,:,j) = RHS * beta_hat;
end

% Choose here the version to save for the baseline results
pred_q = pred_q_V4_AR;
save("../../excel/processed/pred_q.mat", "pred_q"); % time x quantile x country
