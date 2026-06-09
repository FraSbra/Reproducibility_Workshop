function [beta_hat, alpha] = smooth_quantile_reg(y, X, tau)
%QUANTILE_REG Smoothed quantile regression with automatic alpha choice
%
% Inputs:
%   y   - dependent variable (n x 1)
%   X   - regressors (n x k)
%   tau - quantile level, e.g. 0.1, 0.5, 0.9
%
% Outputs:
%   beta_hat - estimated coefficients
%   alpha    - internally selected smoothing parameter

    % Make sure y is a column vector
    y = y(:);

    % Remove invalid rows
    valid = all(isfinite([y X]), 2);
    y = y(valid);
    X = X(valid, :);

    % Initial guess: OLS
    beta0 = (X' * X) \ (X' * y);

    % Initial residuals
    resid0 = y - X * beta0;
    n = length(y);

    % Robust scale estimate based on MAD
    med_resid = median(resid0);
    mad_resid = median(abs(resid0 - med_resid));
    scale_resid = 1.4826 * mad_resid;

    % Fallback if MAD is zero or invalid
    if ~isfinite(scale_resid) || scale_resid <= 0
        scale_resid = std(resid0);
    end

    % Second fallback
    if ~isfinite(scale_resid) || scale_resid <= 0
        scale_resid = std(y);
    end

    % Final fallback
    if ~isfinite(scale_resid) || scale_resid <= 0
        scale_resid = 1;
    end

    % Automatic alpha choice
    % Small enough to approximate quantile loss,
    % but not too small to create numerical problems
    alpha = max(1e-6, 0.1 * scale_resid * n^(-1/3));

    % Smoothed quantile loss
    loss_fun = @(b) smooth_loss(b, y, X, tau, alpha);

    % Optimization options
    options = optimoptions('fminunc', ...
                           'Algorithm', 'quasi-newton', ...
                           'Display', 'off', ...
                           'MaxIterations', 1000, ...
                           'MaxFunctionEvaluations', 5000);

    % Minimize
    beta_hat = fminunc(loss_fun, beta0, options);

end

function f = smooth_loss(b, y, X, tau, alpha)

    u = y - X * b;
    z = -u / alpha;

    % Numerically stable softplus:
    % log(1 + exp(z)) = max(z,0) + log(1 + exp(-abs(z)))
    softplus_z = max(z, 0) + log1p(exp(-abs(z)));

    f = sum(tau * u + alpha * softplus_z);

end