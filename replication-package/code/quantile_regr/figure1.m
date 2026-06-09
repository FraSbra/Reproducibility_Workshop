%% Figure 1: Downside and upside plots

figure;
tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

colors = [0.2 0.4 0.8;   % blue for downside
          0.8 0.2 0.2];  % red for upside

recessions = [datetime(1974,10,1), datetime(1975,3,31);
              datetime(1980,1,1),  datetime(1980,9,30);
              datetime(1982,1,1),  datetime(1982,12,31);
              datetime(1992,1,1),  datetime(1993,3,31);
              datetime(2008,1,1),  datetime(2009,6,30);
              datetime(2011,7,1),  datetime(2013,3,31);
              datetime(2020,1,1),  datetime(2020,6,30)];

for i = 1:size(gdp,2)
    nexttile;
    
    downside = -(pred_q(:,10,i) - pred_q(:,1,i));
    upside   =   pred_q(:,19,i) - pred_q(:,10,i);
    
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
    
    title(countryNames{i}, 'FontSize', 9);
    xlim([datetime(1980,1,1), date(end)]);
    ylim([ymin, ymax]);
    grid on;
    box off;
end

% Hide empty tiles
i = 12;
nexttile;
axis off;


% Single legend at the bottom
lgd = legend([h1, h2], {'Downside uncertainty', 'Upside uncertainty'}, ...
    'Orientation', 'horizontal', 'FontSize', 9);
lgd.Layout.Tile = 'south';

% === SAVE FIGURE AS PNG ===
savePath = fullfile(saveResults_fig, 'down_up_ALL.png');
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, savePath, '-dpng', '-r300');