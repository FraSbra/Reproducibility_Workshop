%% Table 1: In-sample accuracy

pred_q_allVersions = cat(4,pred_q_V1,pred_q_V2);
pred_q_allVersions = cat(4,pred_q_allVersions,pred_q_V3);
pred_q_allVersions = cat(4,pred_q_allVersions,pred_q_V4);
pred_q_allVersions = cat(4,pred_q_allVersions,pred_q_V5);
pred_q_allVersions = cat(4,pred_q_allVersions,pred_q_V6);
pred_q_allVersions = cat(4,pred_q_allVersions,pred_q_V7);
pred_q_allVersions = cat(4,pred_q_allVersions,pred_q_V1_AR);
pred_q_allVersions = cat(4,pred_q_allVersions,pred_q_V2_AR);
pred_q_allVersions = cat(4,pred_q_allVersions,pred_q_V3_AR);
pred_q_allVersions = cat(4,pred_q_allVersions,pred_q_V4_AR);
pred_q_allVersions = cat(4,pred_q_allVersions,pred_q_V5_AR);
pred_q_allVersions = cat(4,pred_q_allVersions,pred_q_V6_AR);
pred_q_allVersions = cat(4,pred_q_allVersions,pred_q_V7_AR);

RMSE = NaN(size(gdp, 2), size(pred_q_allVersions, 4));
MAE = NaN(size(gdp, 2), size(pred_q_allVersions, 4));
QS = NaN(size(gdp, 2), size(pred_q_allVersions, 4));
qwCRPS = NaN(size(gdp, 2), size(pred_q_allVersions, 4));

for j=1:size(pred_q_allVersions,4)
    for i=1:size(gdp,2)
        [RMSE(i,j), MAE(i,j),~] = Point_Forecast_ACCURACY(yy(:,i), pred_q_allVersions(:,:,i,j));
        QS(i,j) = quantile_score(yy(:,i), pred_q_allVersions(:,2,i,j), 0.05);
        qwCRPS(i,j) = weighted_average(yy(:,i), pred_q_allVersions(:,:,i,j), quantiles);
    end
end


results_abs = table(mean(RMSE,1)', mean(MAE,1)', mean(QS,1)', mean(qwCRPS,1)', ...
    'VariableNames', {'RMSE', 'MAE', 'QS', 'qwCRPS'}, ...
    'RowNames',      {'V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V1-AR', 'V2-AR', 'V3-AR', 'V4-AR', 'V5-AR', 'V6-AR', 'V7-AR'});

baseline = results_abs{'V4-AR', :};
results_table = results_abs;
results_table{:, :} = results_abs{:, :} ./ repmat(baseline, height(results_abs), 1);

% --- Store tables as txt ---
fid = fopen(fullfile(saveResults_tab, 'In_Sample_accuracy.txt'), 'w');

% Header
fprintf(fid, '\\begin{table}[htbp]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\begin{tabular}{lcccc}\n');
fprintf(fid, '\\hline\\hline\n');
fprintf(fid, ' & RMSE & MAE & QS & qwCRPS \\\\\n');
fprintf(fid, '\\hline\n');

% Rows
rowNames = results_table.Properties.RowNames;
for i = 1:height(results_table)
    fprintf(fid, '%s & %.4f & %.4f & %.4f & %.4f \\\\\n', ...
        rowNames{i}, ...
        results_table.RMSE(i), ...
        results_table.MAE(i), ...
        results_table.QS(i), ...
        results_table.qwCRPS(i));
end

% Footer
fprintf(fid, '\\hline\\hline\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\caption{In-sample accuracy relative to V4 baseline}\n');
fprintf(fid, '\\label{tab:results}\n');
fprintf(fid, '\\end{table}\n');

fclose(fid);
