function qwCRPS = weighted_average(target, pred_percentiles, percentiles)
    % INPUT: 
    % target: Tx1 target variable of quantile regression
    % pred_percentiles: TxJ grid of predicted percentiles

    % Output:
    % qwCRPS: 1x1 quantile-weighted continous ranked probability score - 
    % - Gneiting & Ranjan (2011)
    Z = [target pred_percentiles]; 
    good = all(~isnan(Z), 2); % Get rid of all nan values
    target = target(good);
    pred_percentiles = pred_percentiles(good,:);
    delta_tau = percentiles(2) - percentiles(1);

    % Choose weight function (tail-focused example)
    w = (2*percentiles - 1).^2;  % emphasize tails
    
    QS_weighted = NaN(length(percentiles), 1);
    for i = 1:length(percentiles)
        u = target - pred_percentiles(:,i);
        QS_weighted(i) = w(i) * mean(u .* (percentiles(i) - (u < 0)));
    end
    
    qwCRPS = 2 * sum(QS_weighted) * delta_tau;





end