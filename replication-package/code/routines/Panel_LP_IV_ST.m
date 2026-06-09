function [Collect_LP]=Panel_LP_IV_ST(YY,XX,Controls,DetermControls,Instr,StateVar,a_id,a_time,settingsLP)
% Panel LP-IV with Smooth Transition (Auerbach-Gorodnichenko style)
%
% Inputs:   YY:             TxN dependent variable (panel outcome)
%           XX:             TxN matrix of regressors OR Tx1 when common across units
%           Controls:       TxNxK 3D matrix of additional controls (panel)
%           DetermControls: TxNxKd 3D matrix of deterministic controls, can be []
%           Instr:          TxN or TxNxQ or TxQ matrix of instruments
%           StateVar:       Tx1 state variable s_t (common across units)
%                           Used to build logistic transition F(s_t)
%           a_id:           TxN matrix of unit identifiers
%           a_time:         TxN matrix of time identifiers
%           settingsLP:     struct with fields:
%               .maxH           = maximum horizon H
%               .MaxLPLagsOwn   = lags of ΔY used as controls
%               .MaxLPLagsOther = lags of other controls
%               .gamma          = slope of logistic function (default: 1.5)
%
% Outputs:  Collect_LP struct:
%           .Coeff_high     = (H+1)x1  IRF coefficients for high-state regime (F(s)~1)
%           .StdErr_high    = (H+1)x1  standard errors, high-state regime
%           .Coeff_low      = (H+1)x1  IRF coefficients for low-state regime (1-F(s)~1)
%           .StdErr_low     = (H+1)x1  standard errors, low-state regime
%           .F_values       = Tx1 transition weights F(s_t) for inspection
%           .gamma          = scalar gamma used
%           .c              = scalar c used (median of StateVar)

%% Transition function parameters

gamma = 1.5;  % default slope
if isfield(settingsLP, 'gamma')
    gamma = settingsLP.gamma;
end

% Standardize state variable (common practice to make gamma scale-free)
s = StateVar(:);                          % Tx1
s_std = (s - nanmean(s)) / nanstd(s);    % standardize

c = nanmedian(s_std);                     % threshold = median (after standardizing, ~0)

% Logistic transition: F(s_t) in (0,1), high values = "high state"
F_vec = 1 ./ (1 + exp(-gamma * (s_std - c)));   % Tx1

%% Settings

maxH          = settingsLP.maxH;
MaxLPLagsOwn  = settingsLP.MaxLPLagsOwn;
MaxLPLagsOther= settingsLP.MaxLPLagsOther;

yy_mat = YY;
[~,NN] = size(yy_mat);

%% Expand common inputs to panel

if size(XX,2) == 1
    XX_mat = repmat(XX,1,NN);
else
    XX_mat = XX;
end

if size(Instr,2) < NN
    Instr_mat = zeros(size(Instr,1), NN, size(Instr,2));
    for ii = 1:size(Instr,2)
        Instr_mat(:,:,ii) = repmat(Instr(:,ii),1,NN);
    end
else
    Instr_mat = Instr;
end

% Lag F_vec by one period: regime must be predetermined w.r.t. the shock
% F(s_{t-1}) avoids simultaneity between the transition weight and dx_t
F_vec_lag = mat_single_leadlag(F_vec, -1);   % Tx1, one lag
F_mat = repmat(F_vec_lag, 1, NN);            % TxN

%% Core variables

vec_id   = a_id(:);
vec_year = a_time(:);

xx   = XX_mat;
dx   = MatrixDiff(xx);
dx   = dx / median(nanstd(dx));   % normalization (same as original)

dy   = yy_mat - mat_single_leadlag(yy_mat,-1);

InstrumentSelect = Instr_mat;

%% Output containers

Collect_LP.Coeff_high  = zeros(maxH+1,1);
Collect_LP.StdErr_high = zeros(maxH+1,1);
Collect_LP.Coeff_low   = zeros(maxH+1,1);
Collect_LP.StdErr_low  = zeros(maxH+1,1);
Collect_LP.F_values    = F_vec_lag;   % F(s_{t-1}), as used in estimation
Collect_LP.gamma       = gamma;
Collect_LP.c           = c;

%% Main LP loop

for hh = 0:maxH

    %% Outcome: cumulative change from t-1 to t+h
    yy_tplus_h     = mat_single_leadlag(yy_mat,hh) - mat_single_leadlag(yy_mat,-1);
    vec_yy_tplus_h = yy_tplus_h(:);

    %% Instruments (vectorized)
    vec_Instr = [];
    for ivij = 1:size(InstrumentSelect,3)
        InstrumentSelect_i = InstrumentSelect(:,:,ivij);
        vec_Instr = [vec_Instr, InstrumentSelect_i(:)];
    end

    %% Regressors: F(s)*dx and (1-F(s))*dx
    vec_dx   = dx(:);
    vec_F    = F_mat(:);          % TN x 1

    vec_Fdx     = vec_F   .* vec_dx;   % high-state weighted shock
    vec_1mFdx   = (1-vec_F) .* vec_dx; % low-state weighted shock

    %% Instruments: F(s)*z and (1-F(s))*z
    vec_Instr_high = vec_F   .* vec_Instr;
    vec_Instr_low  = (1-vec_F) .* vec_Instr;
    vec_Instr_ST   = [vec_Instr_high, vec_Instr_low];

    %% Lag controls
    vec_Controls_t = [];

    for lag = 1:MaxLPLagsOwn
        mat_elem_dy = DemeanVariables(mat_single_leadlag(dy,-lag));
        vec_Controls_t = [vec_Controls_t, mat_elem_dy(:)];
    end

    for pick = 1:size(Controls,3)
        pickControl = DemeanVariables(Controls(:,:,pick));
        for lag = 1:MaxLPLagsOther
            mat_elem_Control = mat_single_leadlag(pickControl,-lag);
            vec_Controls_t = [vec_Controls_t, mat_elem_Control(:)];
        end
    end

    %% Deterministic controls
    if sum(size(DetermControls)) ~= 0
        for pick = 1:size(DetermControls,3)
            if sum(sum(DetermControls(1:end-hh,:,pick))) > 0
                pickDetControl = DetermControls(:,:,pick);
                vec_Controls_t = [vec_Controls_t, pickDetControl(:)];
            else
                fprintf(['Dropped Deterministic Control <<' num2str(pick) '>>, at Horizon [[' num2str(hh) ']]\n']);
            end
        end
    end

    %% IV-FE estimation
    % Endogenous regressors: [F(s)*dx, (1-F(s))*dx] — both instrumented
    % Instruments:           [F(s)*z,  (1-F(s))*z]
    % Controls:              lag controls (treated as exogenous)

    % NaN guard: also catches NaNs propagated from StateVar into F_mat
    Select = isfinite(sum([vec_yy_tplus_h, vec_Fdx, vec_1mFdx, vec_F, vec_Controls_t, vec_Instr_ST], 2));

    ivfe = ivpanel(vec_id(Select,:), vec_year(Select,:), vec_yy_tplus_h(Select,:), ...
                   [vec_Fdx(Select,:), vec_1mFdx(Select,:), vec_Controls_t(Select,:)], ...
                   vec_Instr_ST(Select,:), ...
                   'fe', 'endog', [1, 2]);

    % coef(1) = beta_high (coefficient on F(s)*dx)
    % coef(2) = beta_low  (coefficient on (1-F(s))*dx)
    Collect_LP.Coeff_high(hh+1,1)  = ivfe.coef(1,1);
    Collect_LP.StdErr_high(hh+1,1) = ivfe.stderr(1,1);
    Collect_LP.Coeff_low(hh+1,1)   = ivfe.coef(2,1);
    Collect_LP.StdErr_low(hh+1,1)  = ivfe.stderr(2,1);

    clear ivfe

end