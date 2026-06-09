function [beta_hat, phi, Q_fitted] = smooth_quantile_reg_ar(y, X, tau)
%SMOOTH_QUANTILE_REG_AR Smoothed quantile regression with lagged fitted quantile
%
% Inputs:
%   y   - dependent variable (n x 1)
%   X   - regressors (n x k)
%   tau - quantile level, e.g. 0.1, 0.5, 0.9
%
% Outputs:
%   beta_hat - coefficients on X (k x 1)
%   phi      - coefficient on lagged fitted quantile (scalar)
%   Q_fitted - fitted quantile sequence (n_orig x 1), NaN where y is NaN
%              and at first valid period (no lagged quantile available)
%
% Out-of-sample forecast:
%   Q_next = X_new * beta_hat + phi * Q_fitted(end)

    y = y(:);
    n_orig = length(y);

    % Remove invalid rows (estimation only)
    valid  = all(isfinite([y X]), 2);
    y_est  = y(valid);
    X_est  = X(valid, :);

    n = length(y_est);

    % ----------------------------------------------------------------
    % Step 1: static fit to get Q_lagged
    % ----------------------------------------------------------------
    [beta_static, alpha] = smooth_quantile_reg(y_est, X_est, tau);
    Q_static = X_est * beta_static;         % n x 1 fitted quantiles

    % Lagged fitted quantile (NaN at t=1)
    Q_lagged = [NaN; Q_static(1:end-1)];   % n x 1

    % Drop first row (no lagged quantile available)
    y_aug = y_est(2:end);
    X_aug = [X_est(2:end,:), Q_lagged(2:end)];   % (n-1) x (k+1)

    % ----------------------------------------------------------------
    % Step 2: augmented static fit for warm start
    % ----------------------------------------------------------------
    [theta0, ~] = smooth_quantile_reg(y_aug, X_aug, tau);
    % theta0 is (k+1 x 1): [beta; phi]

    % Seed for the recursion: static fit at t=1
    Q_init = Q_static(1);

    % ----------------------------------------------------------------
    % Step 3: recursive optimization (constrained, |phi| < 1)
    % ----------------------------------------------------------------
    k  = size(X_est, 2);
    lb = [-Inf(k,1); -0.9999];
    ub = [ Inf(k,1);  0.9999];

    loss_fun = @(theta) caviar_smooth_loss(theta, y_est, X_est, tau, alpha, Q_init);

    options = optimoptions('fmincon', ...
                           'Algorithm', 'interior-point', ...
                           'Display', 'off', ...
                           'MaxIterations', 2000, ...
                           'MaxFunctionEvaluations', 10000);

    theta_hat = fmincon(loss_fun, theta0, [], [], [], [], lb, ub, [], options);

    % ----------------------------------------------------------------
    % Step 4: extract coefficients
    % ----------------------------------------------------------------
    beta_hat = theta_hat(1:k);
    phi      = theta_hat(k+1);

    % ----------------------------------------------------------------
    % Step 5: reconstruct fitted quantiles with final parameters
    % ----------------------------------------------------------------
    Q_est    = zeros(n, 1);
    Q_est(1) = Q_init;
    for t = 2:n
        Q_est(t) = X_est(t,:) * beta_hat + phi * Q_est(t-1);
    end

    % Map back to original length, NaN elsewhere
    Q_fitted        = NaN(n_orig, 1);
    Q_fitted(valid) = Q_est;

    % First valid period has no lagged quantile, so set to NaN
    first_valid           = find(valid, 1, 'first');
    Q_fitted(first_valid) = NaN;

end


function f = caviar_smooth_loss(theta, y, X, tau, alpha, Q_init)
%CAVIAR_SMOOTH_LOSS Smoothed quantile loss with recursive quantile propagation

    n = length(y);
    k = size(X, 2);

    beta = theta(1:k);
    phi  = theta(k+1);

    % Recursive propagation of fitted quantiles
    Q    = zeros(n, 1);
    Q(1) = Q_init;
    for t = 2:n
        Q(t) = X(t,:) * beta + phi * Q(t-1);
    end

    % Residuals
    u = y - Q;

    % Smoothed quantile loss
    z          = -u / alpha;
    softplus_z = max(z, 0) + log1p(exp(-abs(z)));
    f          = sum(tau * u + alpha * softplus_z);

end