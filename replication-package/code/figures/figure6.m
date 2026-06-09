%% Figure 6: KBO decomposition with CLIFS -> WORK IN PROGRESS

CLIFS = readtable("data_filtered.xlsx", "Sheet", "CLIFS");
clifs = table2array(CLIFS(:,2:end));
MediatVariable = clifs; % 'Channel' variable

% Shock and instrumented variable
Instrument = ff4;
XX = fedfunds;

interaction = 'No_Interactions';

% Controls
Controls = cat(3, ...                           
    MatrixDiff(repmat(XX,1,No_Countries)), ...      
    MatrixDiff(repmat(100 .* log(other.DEXUSEU),1,No_Countries)), ...
    MatrixDiff(repmat(eonia,1,No_Countries)), ...
    MatrixDiff(repmat(ciss,1,No_Countries)) ...
);

    
[Collect_LP]=Panel_LP_IV_KBO(yy_d,XX,MediatVariable,Controls,[],Instrument,...
    Data.id,Data.quarter,settingsLP,interaction);
% Run the main function

% --- Low RER vulnerability ---
Collect.CoeffMediatDOWN(:) = Collect_LP.CoeffMediatDOWN .* 0.25 / Collect_LP.delta_hat; % IRFs  
Collect.StdErrMediatDOWN(:) = Collect_LP.StdErrMediatDOWN .* 0.25 / Collect_LP.delta_hat; % StdErr

% --- High RER vulnerability --- 
Collect.CoeffMediatUP(:) = Collect_LP.CoeffMediatUP .* 0.25 / Collect_LP.delta_hat; % IRFs  
Collect.StdErrMediatUP(:) = Collect_LP.StdErrMediatUP .* 0.25 / Collect_LP.delta_hat; % StdErr   
    
%% Plot high vs low vulnerability  

h = (0:settingsLP.maxH)';

% --- Extract IRFs and standard errors ---
irf_low  = Collect.CoeffMediatDOWN(:);
se_low   = Collect.StdErrMediatDOWN(:);

irf_high = Collect.CoeffMediatUP(:);
se_high  = Collect.StdErrMediatUP(:);

% --- Confidence bands: LOW vulnerability (red) ---
ub68_low = irf_low + 1.000 * se_low;
lb68_low = irf_low - 1.000 * se_low;
ub90_low = irf_low + 1.645 * se_low;
lb90_low = irf_low - 1.645 * se_low;

% --- Confidence bands: HIGH vulnerability (black) ---
ub68_high = irf_high + 1.000 * se_high;
lb68_high = irf_high - 1.000 * se_high;
ub90_high = irf_high + 1.645 * se_high;
lb90_high = irf_high - 1.645 * se_high;

% --- Colours ---
col_low_90  = [1.00 0.70 0.70];   % light red
col_low_68  = [0.85 0.15 0.15];   % dark  red
col_high_90 = [0.70 0.70 0.70];   % light grey (black family)
col_high_68 = [0.15 0.15 0.15];   % near-black

x_fill = [h; flipud(h)];

figure;
hold on;

% ---- LOW vulnerability (red) ----
fill(x_fill, [ub90_low; flipud(lb90_low)], col_low_90, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.35, 'HandleVisibility', 'off');
fill(x_fill, [ub68_low; flipud(lb68_low)], col_low_68, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.40, 'HandleVisibility', 'off');
p_low = plot(h, irf_low, '-', 'Color', col_low_68, 'LineWidth', 2);

% ---- HIGH vulnerability (black) ----
fill(x_fill, [ub90_high; flipud(lb90_high)], col_high_90, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.35, 'HandleVisibility', 'off');
fill(x_fill, [ub68_high; flipud(lb68_high)], col_high_68, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.40, 'HandleVisibility', 'off');
p_high = plot(h, irf_high, '-', 'Color', col_high_68, 'LineWidth', 2);

% ---- Zero line ----
yline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');

legend([p_low, p_high], ...
       'Low fin. vulnerability', 'High fin. vulnerability', ...
       'Location', 'best');

xlabel('Horizon');
ylabel('Response');
title('IRF to US MP shock – KBO decomposition (CLIFS)');
grid on;
hold off; 
                

%% Plot vulnerability across countries      
MediationRESPONSES = mean(Collect_LP.InteractTermNONdemeaned,1); %
% Average over five horizons of elasticities across countries

[data_sorted, idx] = sort(MediationRESPONSES, 'ascend');
CountryNames = CountryNames(idx);

% Create figure
bb = figure;

% Bar chart: blue bars, no border
b = bar(data_sorted, 'FaceColor', [0 0.4470 0.7410], 'EdgeColor', 'none');

% Set x-axis labels
xticks(1:length(CountryNames));
xticklabels(CountryNames);
xtickangle(40);

% Average line (thick dashed black)
hold on;
avg_val = mean(data_sorted, 'omitnan');
yline(avg_val, '--k', 'LineWidth', 2.5);

% Labels and LaTeX formatting
ylabel('$\%$ deviation from steady state', 'Interpreter', 'latex');

 % Apply LaTeX interpreter to tick labels
ax = gca;
ax.TickLabelInterpreter = 'latex';

% Improve font and spacing
set(gca, 'FontSize', 12);
grid on;

% Set figure size (wider and taller than before)
set(bb, 'Units', 'Inches');
width = 11;  % wider to avoid squeeze
height = 5.5; % more height for clarity
set(bb, 'Position', [1, 1, width, height]);
set(bb, 'PaperPositionMode', 'Auto', ...
    'PaperUnits', 'Inches', ...
    'PaperSize', [width, height]);
