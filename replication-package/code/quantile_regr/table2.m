%% Table 2: Out-of-sample accuracy

start = 150; % Start of out-of-sample testing 

RMSE = NaN(size(gdp,2), size(gdp,1) - start);
MAE = NaN(size(gdp,2), size(gdp,1) - start);
QS = NaN(size(gdp,2), size(gdp,1) - start);
qwCRPS = NaN(size(gdp,2), size(gdp,1) - start);


%% --- V1 ---
for j = 1:size(gdp,2) % Loop over countries
    RHS = [RHS_V1 gdp_lag(:,j)];
    for t = start:(size(pred_q_V1,1)-1) % Loop over out-of-sample observations
        beta_hat = NaN(size(RHS_V1,2) + 1, length(quantiles));
        i = 1;
        for q = quantiles % Loop over quantiles
            RHS_cut = RHS(1:t,:);
            beta_hat(:,i) = smooth_quantile_reg(yy(1:t,j), RHS_cut, q);
            i = i + 1;
        end
        pred_OOS = RHS(t+1,:) * beta_hat; % OOS prediction
        y_true = yy(t+1,j); % 'True' realization
        [~, MAE(j,t-start+1), RMSE(j,t-start+1)] = Point_Forecast_ACCURACY(y_true, pred_OOS);
        QS(j,t-start+1) = quantile_score(y_true, pred_OOS, 0.05);
        qwCRPS(j,t-start+1) = weighted_average(y_true, pred_OOS, quantiles);
    end
end

Z = [RMSE; MAE];
good = all(~isnan(Z), 1);
RMSE = RMSE(:,good);
MAE = MAE(:,good);
QS = QS(:,good);
qwCRPS = qwCRPS(:,good);

RMSE_V1 = sqrt(mean(RMSE(:).^2));
MAE_V1 = mean(MAE(:)); % The same of RMSE by construction
QS_V1 = mean(QS(:));
qwCRPS_V1 = mean(qwCRPS(:));

disp("V1 done");
%% --- V2 ---
for j = 1:size(gdp,2) % Loop over countries
    RHS = [RHS_V2 gdp_lag(:,j)];
    for t = start:(size(pred_q_V2,1)-1) % Loop over out-of-sample observations
        beta_hat = NaN(size(RHS_V2,2) + 1, length(quantiles));
        i = 1;
        for q = quantiles % Loop over quantiles
            RHS_cut = RHS(1:t,:);
            beta_hat(:,i) = smooth_quantile_reg(yy(1:t,j), RHS_cut, q);
            i = i + 1;
        end
        pred_OOS = RHS(t+1,:) * beta_hat; % OOS prediction
        y_true = yy(t+1,j); % 'True' realization
        [~, MAE(j,t-start+1), RMSE(j,t-start+1)] = Point_Forecast_ACCURACY(y_true, pred_OOS);
        QS(j,t-start+1) = quantile_score(y_true, pred_OOS, 0.05);
        qwCRPS(j,t-start+1) = weighted_average(y_true, pred_OOS, quantiles);
    end
end

Z = [RMSE; MAE];
good = all(~isnan(Z), 1);
RMSE = RMSE(:,good);
MAE = MAE(:,good);
QS = QS(:,good);
qwCRPS = qwCRPS(:,good);

RMSE_V2 = sqrt(mean(RMSE(:).^2));
MAE_V2 = mean(MAE(:)); % The same of RMSE by construction
QS_V2 = mean(QS(:));
qwCRPS_V2 = mean(qwCRPS(:));

disp("V2 done");
%% --- V3 ---
for j = 1:size(gdp,2) % Loop over countries
    RHS = [RHS_V3 gdp_lag(:,j)];
    for t = start:(size(pred_q_V3,1)-1) % Loop over out-of-sample observations
        beta_hat = NaN(size(RHS_V3,2) + 1, length(quantiles));
        i = 1;
        for q = quantiles % Loop over quantiles
            RHS_cut = RHS(1:t,:);
            beta_hat(:,i) = smooth_quantile_reg(yy(1:t,j), RHS_cut, q);
            i = i + 1;
        end
        pred_OOS = RHS(t+1,:) * beta_hat; % OOS prediction
        y_true = yy(t+1,j); % 'True' realization
        [~, MAE(j,t-start+1), RMSE(j,t-start+1)] = Point_Forecast_ACCURACY(y_true, pred_OOS);
        QS(j,t-start+1) = quantile_score(y_true, pred_OOS, 0.05);
        qwCRPS(j,t-start+1) = weighted_average(y_true, pred_OOS, quantiles);
    end
end

Z = [RMSE; MAE];
good = all(~isnan(Z), 1);
RMSE = RMSE(:,good);
MAE = MAE(:,good);
QS = QS(:,good);
qwCRPS = qwCRPS(:,good);

RMSE_V3 = sqrt(mean(RMSE(:).^2));
MAE_V3 = mean(MAE(:)); % The same of RMSE by construction
QS_V3 = mean(QS(:));
qwCRPS_V3 = mean(qwCRPS(:));

disp("V3 done");
%% --- V4 ---
for j = 1:size(gdp,2) % Loop over countries
    RHS = [RHS_V4 gdp_lag(:,j)];
    for t = start:(size(pred_q_V4,1)-1) % Loop over out-of-sample observations
        beta_hat = NaN(size(RHS_V4,2) + 1, length(quantiles));
        i = 1;
        for q = quantiles % Loop over quantiles
            RHS_cut = RHS(1:t,:);
            beta_hat(:,i) = smooth_quantile_reg(yy(1:t,j), RHS_cut, q);
            i = i + 1;
        end
        pred_OOS = RHS(t+1,:) * beta_hat; % OOS prediction
        y_true = yy(t+1,j); % 'True' realization
        [~, MAE(j,t-start+1), RMSE(j,t-start+1)] = Point_Forecast_ACCURACY(y_true, pred_OOS);
        QS(j,t-start+1) = quantile_score(y_true, pred_OOS, 0.05);
        qwCRPS(j,t-start+1) = weighted_average(y_true, pred_OOS, quantiles);
    end
end

Z = [RMSE; MAE];
good = all(~isnan(Z), 1);
RMSE = RMSE(:,good);
MAE = MAE(:,good);
QS = QS(:,good);
qwCRPS = qwCRPS(:,good);

RMSE_V4 = sqrt(mean(RMSE(:).^2));
MAE_V4 = mean(MAE(:)); % The same of RMSE by construction
QS_V4 = mean(QS(:));
qwCRPS_V4 = mean(qwCRPS(:));

disp("V4 done");
%% --- V5 ---
for j = 1:size(gdp,2) % Loop over countries
    RHS = [RHS_V2 debt_lag(:,j) gdp_lag(:,j)];
    for t = start:(size(pred_q_V5,1)-1) % Loop over out-of-sample observations
        beta_hat = NaN(size(RHS_V2,2) + 2, length(quantiles));
        i = 1;
        for q = quantiles % Loop over quantiles
            RHS_cut = RHS(1:t,:);
            beta_hat(:,i) = smooth_quantile_reg(yy(1:t,j), RHS_cut, q);
            i = i + 1;
        end
        pred_OOS = RHS(t+1,:) * beta_hat; % OOS prediction
        y_true = yy(t+1,j); % 'True' realization
        [~, MAE(j,t-start+1), RMSE(j,t-start+1)] = Point_Forecast_ACCURACY(y_true, pred_OOS);
        QS(j,t-start+1) = quantile_score(y_true, pred_OOS, 0.05);
        qwCRPS(j,t-start+1) = weighted_average(y_true, pred_OOS, quantiles);
    end
end

Z = [RMSE; MAE];
good = all(~isnan(Z), 1);
RMSE = RMSE(:,good);
MAE = MAE(:,good);
QS = QS(:,good);
qwCRPS = qwCRPS(:,good);

RMSE_V5 = sqrt(mean(RMSE(:).^2));
MAE_V5 = mean(MAE(:)); % The same of RMSE by construction
QS_V5 = mean(QS(:));
qwCRPS_V5 = mean(qwCRPS(:));

disp("V5 done");
%% --- V6 ---
for j = 1:size(gdp,2) % Loop over countries
    RHS = [RHS_V6 gdp_lag(:,j)];
    for t = start:(size(pred_q_V6,1)-1) % Loop over out-of-sample observations
        beta_hat = NaN(size(RHS_V6,2) + 1, length(quantiles));
        i = 1;
        for q = quantiles % Loop over quantiles
            RHS_cut = RHS(1:t,:);
            beta_hat(:,i) = smooth_quantile_reg(yy(1:t,j), RHS_cut, q);
            i = i + 1;
        end
        pred_OOS = RHS(t+1,:) * beta_hat; % OOS prediction
        y_true = yy(t+1,j); % 'True' realization
        [~, MAE(j,t-start+1), RMSE(j,t-start+1)] = Point_Forecast_ACCURACY(y_true, pred_OOS);
        QS(j,t-start+1) = quantile_score(y_true, pred_OOS, 0.05);
        qwCRPS(j,t-start+1) = weighted_average(y_true, pred_OOS, quantiles);
    end
end

Z = [RMSE; MAE];
good = all(~isnan(Z), 1);
RMSE = RMSE(:,good);
MAE = MAE(:,good);
QS = QS(:,good);
qwCRPS = qwCRPS(:,good);

RMSE_V6 = sqrt(mean(RMSE(:).^2));
MAE_V6 = mean(MAE(:)); % The same of RMSE by construction
QS_V6 = mean(QS(:));
qwCRPS_V6 = mean(qwCRPS(:));

disp("V6 done");
%% --- V7 ---
for j = 1:size(gdp,2) % Loop over countries
    RHS = [RHS_V7 gdp_lag(:,j)];
    for t = start:(size(pred_q_V7,1)-1) % Loop over out-of-sample observations
        beta_hat = NaN(size(RHS_V7,2) + 1, length(quantiles));
        i = 1;
        for q = quantiles % Loop over quantiles
            RHS_cut = RHS(1:t,:);
            beta_hat(:,i) = smooth_quantile_reg(yy(1:t,j), RHS_cut, q);
            i = i + 1;
        end
        pred_OOS = RHS(t+1,:) * beta_hat; % OOS prediction
        y_true = yy(t+1,j); % 'True' realization
        [~, MAE(j,t-start+1), RMSE(j,t-start+1)] = Point_Forecast_ACCURACY(y_true, pred_OOS);
        QS(j,t-start+1) = quantile_score(y_true, pred_OOS, 0.05);
        qwCRPS(j,t-start+1) = weighted_average(y_true, pred_OOS, quantiles);
    end
end

Z = [RMSE; MAE];
good = all(~isnan(Z), 1);
RMSE = RMSE(:,good);
MAE = MAE(:,good);
QS = QS(:,good);
qwCRPS = qwCRPS(:,good);

RMSE_V7 = sqrt(mean(RMSE(:).^2));
MAE_V7 = mean(MAE(:)); % The same of RMSE by construction
QS_V7 = mean(QS(:));
qwCRPS_V7 = mean(qwCRPS(:));

disp("V7 done");
%% --- V1-AR ---
for j = 1:size(gdp,2) % Loop over countries
    RHS = [RHS_V1 gdp_lag(:,j)];
    for t = start:(size(pred_q_V1,1)-1) % Loop over out-of-sample observations
        beta_hat = NaN(size(RHS_V1,2) + 1, length(quantiles));
        phi = NaN(1,length(quantiles));
        i = 1;
        for q = quantiles % Loop over quantiles
            RHS_cut = RHS(1:t,:);
            [beta_hat(:,i), phi(i), pred_q_V1_AR(1:t,i,j)] = smooth_quantile_reg_ar(yy(1:t,j), RHS_cut, q);
            i = i + 1;
        end
        pred_OOS =  RHS(t+1,:) * beta_hat + phi .* pred_q_V1_AR(t,:,j); % OOS prediction
        y_true = yy(t+1,j); % 'True' realization
        [~, MAE(j,t-start+1), RMSE(j,t-start+1)] = Point_Forecast_ACCURACY(y_true, pred_OOS);
        QS(j,t-start+1) = quantile_score(y_true, pred_OOS, 0.05);
        qwCRPS(j,t-start+1) = weighted_average(y_true, pred_OOS, quantiles);
    end
end

Z = [RMSE; MAE];
good = all(~isnan(Z), 1);
RMSE = RMSE(:,good);
MAE = MAE(:,good);
QS = QS(:,good);
qwCRPS = qwCRPS(:,good);

RMSE_V1_AR = sqrt(mean(RMSE(:).^2));
MAE_V1_AR = mean(MAE(:)); % The same of RMSE by construction
QS_V1_AR = mean(QS(:));
qwCRPS_V1_AR = mean(qwCRPS(:));

disp("V1-AR done");
%% --- V2-AR ---
for j = 1:size(gdp,2) % Loop over countries
    RHS = [RHS_V2 gdp_lag(:,j)];
    for t = start:(size(pred_q_V2,1)-1) % Loop over out-of-sample observations
        beta_hat = NaN(size(RHS_V2,2) + 1, length(quantiles));
        phi = NaN(1,length(quantiles));
        i = 1;
        for q = quantiles % Loop over quantiles
            RHS_cut = RHS(1:t,:);
            [beta_hat(:,i), phi(i), pred_q_V2_AR(1:t,i,j)] = smooth_quantile_reg_ar(yy(1:t,j), RHS_cut, q);
            i = i + 1;
        end
        pred_OOS =  RHS(t+1,:) * beta_hat + phi .* pred_q_V2_AR(t,:,j); % OOS prediction
        y_true = yy(t+1,j); % 'True' realization
        [~, MAE(j,t-start+1), RMSE(j,t-start+1)] = Point_Forecast_ACCURACY(y_true, pred_OOS);
        QS(j,t-start+1) = quantile_score(y_true, pred_OOS, 0.05);
        qwCRPS(j,t-start+1) = weighted_average(y_true, pred_OOS, quantiles);
    end
end

Z = [RMSE; MAE];
good = all(~isnan(Z), 1);
RMSE = RMSE(:,good);
MAE = MAE(:,good);
QS = QS(:,good);
qwCRPS = qwCRPS(:,good);

RMSE_V2_AR = sqrt(mean(RMSE(:).^2));
MAE_V2_AR = mean(MAE(:)); % The same of RMSE by construction
QS_V2_AR = mean(QS(:));
qwCRPS_V2_AR = mean(qwCRPS(:));

disp("V2-AR done");
%% --- V3-AR ---
for j = 1:size(gdp,2) % Loop over countries
    RHS = [RHS_V3 gdp_lag(:,j)];
    for t = start:(size(pred_q_V3,1)-1) % Loop over out-of-sample observations
        beta_hat = NaN(size(RHS_V3,2) + 1, length(quantiles));
        phi = NaN(1,length(quantiles));
        i = 1;
        for q = quantiles % Loop over quantiles
            RHS_cut = RHS(1:t,:);
            [beta_hat(:,i), phi(i), pred_q_V3_AR(1:t,i,j)] = smooth_quantile_reg_ar(yy(1:t,j), RHS_cut, q);
            i = i + 1;
        end
        pred_OOS =  RHS(t+1,:) * beta_hat + phi .* pred_q_V3_AR(t,:,j); % OOS prediction
        y_true = yy(t+1,j); % 'True' realization
        [~, MAE(j,t-start+1), RMSE(j,t-start+1)] = Point_Forecast_ACCURACY(y_true, pred_OOS);
        QS(j,t-start+1) = quantile_score(y_true, pred_OOS, 0.05);
        qwCRPS(j,t-start+1) = weighted_average(y_true, pred_OOS, quantiles);
    end
end

Z = [RMSE; MAE];
good = all(~isnan(Z), 1);
RMSE = RMSE(:,good);
MAE = MAE(:,good);
QS = QS(:,good);
qwCRPS = qwCRPS(:,good);

RMSE_V3_AR = sqrt(mean(RMSE(:).^2));
MAE_V3_AR = mean(MAE(:)); % The same of RMSE by construction
QS_V3_AR = mean(QS(:));
qwCRPS_V3_AR = mean(qwCRPS(:));

disp("V3-AR done");
%% --- V4-AR ---
for j = 1:size(gdp,2) % Loop over countries
    RHS = [RHS_V4 gdp_lag(:,j)];
    for t = start:(size(pred_q_V4,1)-1) % Loop over out-of-sample observations
        beta_hat = NaN(size(RHS_V4,2) + 1, length(quantiles));
        phi = NaN(1,length(quantiles));
        i = 1;
        for q = quantiles % Loop over quantiles
            RHS_cut = RHS(1:t,:);
            [beta_hat(:,i), phi(i), pred_q_V4_AR(1:t,i,j)] = smooth_quantile_reg_ar(yy(1:t,j), RHS_cut, q);
            i = i + 1;
        end
        pred_OOS =  RHS(t+1,:) * beta_hat + phi .* pred_q_V4_AR(t,:,j); % OOS prediction
        y_true = yy(t+1,j); % 'True' realization
        [~, MAE(j,t-start+1), RMSE(j,t-start+1)] = Point_Forecast_ACCURACY(y_true, pred_OOS);
        QS(j,t-start+1) = quantile_score(y_true, pred_OOS, 0.05);
        qwCRPS(j,t-start+1) = weighted_average(y_true, pred_OOS, quantiles);
    end
end

Z = [RMSE; MAE];
good = all(~isnan(Z), 1);
RMSE = RMSE(:,good);
MAE = MAE(:,good);
QS = QS(:,good);
qwCRPS = qwCRPS(:,good);

RMSE_V4_AR = sqrt(mean(RMSE(:).^2));
MAE_V4_AR = mean(MAE(:)); % The same of RMSE by construction
QS_V4_AR = mean(QS(:));
qwCRPS_V4_AR = mean(qwCRPS(:));

disp("V4-AR done");
%% --- V5-AR ---
for j = 1:size(gdp,2) % Loop over countries
    RHS = [RHS_V2 debt_lag(:,j) gdp_lag(:,j)];
    for t = start:(size(pred_q_V5,1)-1) % Loop over out-of-sample observations
        beta_hat = NaN(size(RHS_V2,2) + 2, length(quantiles));
        phi = NaN(1,length(quantiles));
        i = 1;
        for q = quantiles % Loop over quantiles
            RHS_cut = RHS(1:t,:);
            [beta_hat(:,i), phi(i), pred_q_V5_AR(1:t,i,j)] = smooth_quantile_reg_ar(yy(1:t,j), RHS_cut, q);
            i = i + 1;
        end
        pred_OOS =  RHS(t+1,:) * beta_hat + phi .* pred_q_V5_AR(t,:,j); % OOS prediction
        y_true = yy(t+1,j); % 'True' realization
        [~, MAE(j,t-start+1), RMSE(j,t-start+1)] = Point_Forecast_ACCURACY(y_true, pred_OOS);
        QS(j,t-start+1) = quantile_score(y_true, pred_OOS, 0.05);
        qwCRPS(j,t-start+1) = weighted_average(y_true, pred_OOS, quantiles);
    end
end

Z = [RMSE; MAE];
good = all(~isnan(Z), 1);
RMSE = RMSE(:,good);
MAE = MAE(:,good);
QS = QS(:,good);
qwCRPS = qwCRPS(:,good);

RMSE_V5_AR = sqrt(mean(RMSE(:).^2));
MAE_V5_AR = mean(MAE(:)); % The same of RMSE by construction
QS_V5_AR = mean(QS(:));
qwCRPS_V5_AR = mean(qwCRPS(:));

disp("V5-AR done");
%% --- V6-AR ---
for j = 1:size(gdp,2) % Loop over countries
    RHS = [RHS_V6 gdp_lag(:,j)];
    for t = start:(size(pred_q_V6,1)-1) % Loop over out-of-sample observations
        beta_hat = NaN(size(RHS_V6,2) + 1, length(quantiles));
        phi = NaN(1,length(quantiles));
        i = 1;
        for q = quantiles % Loop over quantiles
            RHS_cut = RHS(1:t,:);
            [beta_hat(:,i), phi(i), pred_q_V6_AR(1:t,i,j)] = smooth_quantile_reg_ar(yy(1:t,j), RHS_cut, q);
            i = i + 1;
        end
        pred_OOS =  RHS(t+1,:) * beta_hat + phi .* pred_q_V6_AR(t,:,j); % OOS prediction
        y_true = yy(t+1,j); % 'True' realization
        [~, MAE(j,t-start+1), RMSE(j,t-start+1)] = Point_Forecast_ACCURACY(y_true, pred_OOS);
        QS(j,t-start+1) = quantile_score(y_true, pred_OOS, 0.05);
        qwCRPS(j,t-start+1) = weighted_average(y_true, pred_OOS, quantiles);
    end
end

Z = [RMSE; MAE];
good = all(~isnan(Z), 1);
RMSE = RMSE(:,good);
MAE = MAE(:,good);
QS = QS(:,good);
qwCRPS = qwCRPS(:,good);

RMSE_V6_AR = sqrt(mean(RMSE(:).^2));
MAE_V6_AR = mean(MAE(:)); % The same of RMSE by construction
QS_V6_AR = mean(QS(:));
qwCRPS_V6_AR = mean(qwCRPS(:));

disp("V6-AR done");
%% --- V7-AR ---
for j = 1:size(gdp,2) % Loop over countries
    RHS = [RHS_V7 gdp_lag(:,j)];
    for t = start:(size(pred_q_V7,1)-1) % Loop over out-of-sample observations
        beta_hat = NaN(size(RHS_V7,2) + 1, length(quantiles));
        phi = NaN(1,length(quantiles));
        i = 1;
        for q = quantiles % Loop over quantiles
            RHS_cut = RHS(1:t,:);
            [beta_hat(:,i), phi(i), pred_q_V7_AR(1:t,i,j)] = smooth_quantile_reg_ar(yy(1:t,j), RHS_cut, q);
            i = i + 1;
        end
        pred_OOS =  RHS(t+1,:) * beta_hat + phi .* pred_q_V7_AR(t,:,j); % OOS prediction
        y_true = yy(t+1,j); % 'True' realization
        [~, MAE(j,t-start+1), RMSE(j,t-start+1)] = Point_Forecast_ACCURACY(y_true, pred_OOS);
        QS(j,t-start+1) = quantile_score(y_true, pred_OOS, 0.05);
        qwCRPS(j,t-start+1) = weighted_average(y_true, pred_OOS, quantiles);
    end
end

Z = [RMSE; MAE];
good = all(~isnan(Z), 1);
RMSE = RMSE(:,good);
MAE = MAE(:,good);
QS = QS(:,good);
qwCRPS = qwCRPS(:,good);

RMSE_V7_AR = sqrt(mean(RMSE(:).^2));
MAE_V7_AR = mean(MAE(:)); % The same of RMSE by construction
QS_V7_AR = mean(QS(:));
qwCRPS_V7_AR = mean(qwCRPS(:));

disp("V7-AR done");
%% Final Table
results_table_OOS_abs = table( ...
    [RMSE_V1;    RMSE_V2;    RMSE_V3;    RMSE_V4;    RMSE_V5;    RMSE_V6;    RMSE_V7;    RMSE_V1_AR;    RMSE_V2_AR;    RMSE_V3_AR;    RMSE_V4_AR;    RMSE_V5_AR;    RMSE_V6_AR;    RMSE_V7_AR], ...
    [MAE_V1;     MAE_V2;     MAE_V3;     MAE_V4;     MAE_V5;     MAE_V6;     MAE_V7;     MAE_V1_AR;     MAE_V2_AR;     MAE_V3_AR;     MAE_V4_AR;     MAE_V5_AR;     MAE_V6_AR;     MAE_V7_AR], ...
    [QS_V1;      QS_V2;      QS_V3;      QS_V4;      QS_V5;      QS_V6;      QS_V7;      QS_V1_AR;      QS_V2_AR;      QS_V3_AR;      QS_V4_AR;      QS_V5_AR;      QS_V6_AR;      QS_V7_AR], ...
    [qwCRPS_V1;  qwCRPS_V2;  qwCRPS_V3;  qwCRPS_V4;  qwCRPS_V5;  qwCRPS_V6;  qwCRPS_V7;  qwCRPS_V1_AR;  qwCRPS_V2_AR;  qwCRPS_V3_AR;  qwCRPS_V4_AR;  qwCRPS_V5_AR;  qwCRPS_V6_AR;  qwCRPS_V7_AR], ...
    'VariableNames', {'RMSE', 'MAE', 'QS', 'qwCRPS'}, ...
    'RowNames',      {'V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V1-AR', 'V2-AR', 'V3-AR', 'V4-AR', 'V5-AR', 'V6-AR', 'V7-AR'});

baseline = results_table_OOS_abs{'V4-AR', :};
results_table_OOS = results_table_OOS_abs;
results_table_OOS{:, :} = results_table_OOS_abs{:, :} ./ repmat(baseline, height(results_table_OOS_abs), 1);

fid = fopen(fullfile(saveResults_tab, 'OOS_accuracy.txt'), 'w');
% Header
fprintf(fid, '\\begin{table}[htbp]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\begin{tabular}{lcccc}\n');
fprintf(fid, '\\hline\\hline\n');
fprintf(fid, ' & RMSE & MAE & QS & qwCRPS \\\\\n');
fprintf(fid, '\\hline\n');
% Rows
rowNames = results_table_OOS.Properties.RowNames;
for i = 1:height(results_table_OOS)
    fprintf(fid, '%s & %.4f & %.4f & %.4f & %.4f \\\\\n', ...
        rowNames{i}, ...
        results_table_OOS.RMSE(i), ...
        results_table_OOS.MAE(i), ...
        results_table_OOS.QS(i), ...
        results_table_OOS.qwCRPS(i));
end
% Footer
fprintf(fid, '\\hline\\hline\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\caption{Out-of-sample accuracy relative to V4-AR baseline}\n');
fprintf(fid, '\\label{tab:oos}\n');
fprintf(fid, '\\end{table}\n');
fclose(fid);
