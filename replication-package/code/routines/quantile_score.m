function QS = quantile_score(target, pred_quantile, tau)
    % INPUT: 
    % target: Tx1 target variable of quantile regression
    % pred_quantile: Tx1 predicted quantile

    % OUTPUT:
    % QS: 1x1 quantile score
    Z = [target pred_quantile]; 
    good = all(~isnan(Z), 2); % Get rid of all nan values
    target = target(good);
    pred_quantile = pred_quantile(good);

    u = target - pred_quantile;
    QS = mean(u .* (tau - (u < 0)));

end
