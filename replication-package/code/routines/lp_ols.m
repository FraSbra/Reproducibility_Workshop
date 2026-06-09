function res = lp_ols(y, X, S, H, hStart, lpMode, seType, c)
% Inputs:   y: Tx1 dependent variable
%           X: Txk matrix of regressors 
%           S: TxM matrix of shocks
%           hStart: integer where to start the LP (usually 0 or 1)
%           lpMode: "cum" --> cumulative changes (i.e. long differences)
%                   "level" --> levels
%           seType: "ols" --> OLS standard errors ('classic')
%                   "hac" --> Newey-West standard errors
%           c: integer 0 --> no constant added
%                      1 --> constant added 
%
% Outputs:  res = struct()
%           res.beta = (H+1-hStart)x1 coefficients of shocks
%           res.se = (H+1-hStart)x1 standard errors of coefficients of shocks
%           res.h = (H+1-hStart)x1 horizons of IRFs

    [T, M] = size(S); % Store dimensions of matrix of shocks 

    beta = NaN(H+1-hStart, M); % Initialize matrix to store coefficients of
    % of shocks (i.e. IRFs)
    se   = NaN(H+1-hStart, M); % Same for corresponding standard errors

    for h = hStart:H % Loop over horizons

        switch lpMode % Build dependent variable for both specifications

            case "level" % Level specification
                y_lead = y(h+1:end);
                yh = [y_lead; NaN(h,1)];

            case "cum" % Long diff specification 
                y_lead = [y(1+h:end); NaN(h,1)];       % y_{t+h}
                y_lag  = [NaN(1,1); y(1:end-1)];       % y_{t-1}
                yh = y_lead - y_lag;

            otherwise
                error("lpMode must be 'level' or 'cum'."); % Sanity check 
        end

        Z = [yh X S]; 
        good = all(~isnan(Z), 2); % Get rid of all nan values

        yreg = yh(good); % Final dependent variable

        if c == 1 % Add constant term in this case
            Xreg = [ones(size(X(good,:),1),1) X(good,:)]; % Final matrix of 
            % controls (case constant)
        elseif c == 0 % Do not add constant term in this case
            Xreg = X(good,:); % Final matrix of controls (case no constant)
        else
            error('c must be 0 or 1') % Sanity check
        end

        Sreg = S(good, :); % Final matrix of shocks

        RHS = [Xreg Sreg]; % Right hand side of LP regression (controls + shocks)

        row = h + 1 - hStart; % Current element of the IRF

        switch string(seType) % Estimate coefficients and standard errors 
            % for 'classic' OLS and Newey-West

            case "ols" % 'Classic' OLS case
                [bhat, ~, r, ~, ~] = regress(yreg, RHS);   % r = residuals

                beta_h = bhat(end-M+1:end);
                beta(row, :) = beta_h';

                XtXinv = inv(RHS' * RHS);

                uhat = r;                            % residuals
                meat = RHS' * ( (uhat.^2) .* RHS );  % X' * diag(uhat.^2) * X
                V_white = XtXinv * meat * XtXinv;    % HC0 / White correction for 
                % heteroskedasticity

                se(row, :) = sqrt(diag(V_white(end-M+1:end, end-M+1:end))).';

            case "hac" % Newey-West case
                L = h; % Number of lags for which we control autocorrelation
                % with this choice, we set it equal to the number of
                % horizons

                if c == 1 % hac function already includes contant --> to 
                    % avoid multicollinearity we have to eliminate it if we
                    % added before
                    RHS = RHS(:,2:end); % Delete the column of the constant
                end

                [~, se_all, coeff] = hac(RHS, yreg, ...
                                         Bandwidth=L+1, ...
                                         Display="off");

                beta_h = coeff(end-M+1:end);
                beta(row, :) = beta_h'; % Store coeffcients of interest (shocks)

                se(row, :) = se_all(end-M+1:end)'; % Store standard errors of
                % coeffcients of interest (shocks)

            otherwise
                error("seType must be 'ols' or 'hac'."); % Sanity check
        end
    end

    res = struct(); % Store results
    res.beta = beta; % IRFs
    res.se   = se; % Standard errors
    res.h    = (hStart:H)'; % Horizons of IRFs

end