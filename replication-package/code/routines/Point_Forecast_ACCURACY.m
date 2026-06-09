function [RMSE,MAE,error] = Point_Forecast_ACCURACY(target,pred_percentiles)
    % INPUT: 
    % target: Tx1 target variable of quantile regression
    % pred_percentiles: TxJ grid of predicted percentiles

    % Output:
    % RMSE: 1x1 Root Mean Squared Error, sqrt( (1/T) * sum_{t=1}^{T} (y_t - y_hat_t)^2 )
    % MAE: 1x1 Mean Absolute Error, (1/T) * sum_{t=1}^{T} |y_t - y_hat_t|

    Z = [target pred_percentiles]; 
    good = all(~isnan(Z), 2); % Get rid of all nan values
    target = target(good);
    pred_percentiles = pred_percentiles(good,:);

    J = size(pred_percentiles, 2);
    
    y_hat = (1/J) .* sum(pred_percentiles, 2); % Fitted mean
  
    error = (1/length(target)) * sum(target - y_hat);
    RMSE = sqrt((1/length(target)) * sum((target - y_hat).^2));
    MAE = (1/length(target)) * sum(abs(target - y_hat));
end