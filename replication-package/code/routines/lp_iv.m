function res = lp_iv(y, X, S, Ziv, H, hStart, lpMode, c)
% Inputs:   y: Tx1 dependent variable
%           X: Txk matrix of regressors 
%           S: Tx1 matrix of endogenous shocks (to be instrumented)
%           Ziv: TxQ matrix of instruments for S (e.g., RR shocks)
%           H: integer maximum number of horizons
%           hStart: integer where to start the LP (usually 0 or 1)
%           lpMode: "cum" --> cumulative changes (i.e. long differences)
%                   "lev" --> levels
%           c: integer 0 --> no constant added
%                      1 --> constant added 
%
% Outputs:  res = struct()
%           res.beta = (H+1-hStart)x1 coefficients of shocks
%           res.se = (H+1-hStart)x1 standard errors of coefficients of shocks
%           res.h = (H+1-hStart)x1 horizons of IRFs

    [T, M] = size(S); % Store dimensions of matrix of shocks  

    beta = NaN(H+1-hStart, M); % Initialize matrix to store coefficients of
    % endogenous shocks (i.e. IV-IRFs)
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

        Zall = [yh X S Ziv]; % Stack all variables needed at horizon h
        good = all(~isnan(Zall), 2); % Get rid of all nan values

        yreg = yh(good); % Final dependent variable

        if c == 1 % Add constant term in this case
            Xreg = [ones(size(X(good,:),1),1) X(good,:)]; % Final matrix of
            % controls (case constant)
        elseif c == 0 % Do not add constant term in this case
            Xreg = X(good,:); % Final matrix of controls (case no constant)
        else
            error('c must be 0 or 1') % Sanity check
        end

        Sreg   = S(good);   % Final matrix of endogenous shocks
        Zivreg = Ziv(good, :); % Final matrix of instruments

        row = h + 1 - hStart; % Current element of the IRF

       [betaiv, seiv] = ivreg(yreg, Sreg, Zivreg, Xreg); % See Michal Kolesár 
       % notes for details

       beta(row, :) = betaiv(2); % Two-stages least squares estimator
       se(row, :) =  seiv(2,2); % Heteroscedasticity-robust standard errors
    end

    res = struct(); % Store results
    res.beta = beta; % IRFs
    res.se   = se; % Standard errors
    res.h    = (hStart:H)'; % Horizons of IRFs

end