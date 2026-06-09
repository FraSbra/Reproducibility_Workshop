function [Collect_LP]=Panel_LP_IV_KBO(YY,XX,MediatVariable,Controls,DetermControls, ...
    Instr,a_id,a_time,settingsLP, interaction)
% Inputs:   YY: TxN dependent variable for the final LP regression (panel outcome)
%           XX: TxN matrix of regressors OR Tx1 when common across units
%               (interpreted as the shock/intervention variable whose first-difference
%               is instrumented in the second-stage LP)
%           MediatVariable: TxN panel variable used to estimate cross-sectional
%               heterogeneity in the shock effect across units
%           Controls: TxNxK 3D matrix of additional controls (panel)
%           DetermControls: TxNxKd 3D matrix of deterministic controls (panel),
%                          can be [] if none
%           Instr:  TxN matrix of instruments OR TxNxQ 3D matrix of instruments
%                  OR TxQ when common across units
%           a_id:   TxN matrix of unit identifiers (panel structure)
%           a_time: TxN matrix of time identifiers (panel structure)
%           settingsLP: struct with LP settings:
%               settingsLP.maxH           = maximum horizon H (integer)
%               settingsLP.MaxLPLagsOwn   = number of lags of DeltaY used as controls (integer)
%               settingsLP.MaxLPLagsOther = number of lags of other controls (integer)
%
% Outputs:  Collect_LP = struct()
%
%           % Step 1: estimated unit-level interaction terms
%           Collect_LP.InteractTerm            = (H+1)xN demeaned unit-specific
%                                               interaction coefficients
%           Collect_LP.InteractTermNONdemeaned = (H+1)xN raw unit-specific
%                                               interaction coefficients
%           Collect_LP.CollectMeanunits        = (H+1)x1 cross-sectional mean of
%                                               unit-specific coefficients
%           Collect_LP.CollectSDTunits         = (H+1)x1 cross-sectional std. dev.
%                                               of unit-specific coefficients
%
%           % Step 2: LP-IV with mediator/interacted shock
%           Collect_LP.ATE            = (H+1)x1 baseline average treatment effect
%                                       of the shock
%           Collect_LP.ATEStdErr      = (H+1)x1 standard errors of baseline effect
%
%           Collect_LP.Mediat         = (H+1)x1 coefficient on the mediator interaction
%                                       (shock x estimated heterogeneity term)
%           Collect_LP.MediatStdErr   = (H+1)x1 standard errors of mediator coefficient
%           Collect_LP.MediatSignificant = (H+1)x1 indicator equal to 1 when the
%                                       mediator coefficient is significant at the
%                                       chosen significance level
%
%           Collect_LP.CoeffMediatDOWN = (H+1)x1 IRF evaluated at a low value of the
%                                        mediator (-1.5 cross-sectional std. dev.)
%           Collect_LP.StdErrMediatDOWN= (H+1)x1 standard error for DOWN IRF
%
%           Collect_LP.CoeffMediatUP   = (H+1)x1 IRF evaluated at a high value of the
%                                        mediator (+1.5 cross-sectional std. dev.)
%           Collect_LP.StdErrMediatUP  = (H+1)x1 standard error for UP IRF

%% SPECIFY TYPE 

KBOType.CaseInteraction = interaction; 
NoInteraction = size(Controls,3)+1;
[T,No_Countries] = size(YY); 

%% If XX Common Across Countries -- Create Panel 
if size(XX,2) == 1
    XX_mat = repmat(XX,1,No_Countries); 
else 
    XX_mat = XX; 
end

%% If Instr Common Across Countries -- Create Panel 
if size(Instr,2) < No_Countries
    for ii = 1:size(Instr,2)
        Instr_mat(:,:,ii) = repmat(Instr(:,ii),1,No_Countries);
    end
else
    Instr_mat = Instr; 
end

%% ESTIMTE THE INTERACTING EFFECTS 

% A) ORGANIZE DATA

yy = MediatVariable; 
dy = yy - mat_single_leadlag(yy,-1); 

xx = XX_mat; 
dx = xx - mat_single_leadlag(xx,-1); 

SignShock = 1;
dx = SignShock*dx./repmat(nanstd(dx),size(dx,1),1);

InstrumentSelect = Instr_mat; 

Countries_f_it = zeros(T*No_Countries,No_Countries);
Countries_Instr_f_it = zeros(T*No_Countries,No_Countries,size(InstrumentSelect,3));

for ii = 1:No_Countries 
    dummy = zeros(T,No_Countries); dummy(:,ii) = dx(:,ii);
    Countries_f_it(:,ii) = dummy(:);
    for jjii = 1:size(Countries_Instr_f_it,3)
    dummy = zeros(T,No_Countries); dummy(:,ii) = squeeze(InstrumentSelect(:,ii,jjii));
    Countries_Instr_f_it(:,ii,jjii) = dummy(:);
    end
end

% B) ESTIMATE PANEL LP 
vec_id = a_id(:);  
vec_year = a_time(:); 

InteractCoeff = NaN*zeros(settingsLP.maxH+1,No_Countries);
InteractCoeffSE = NaN*zeros(settingsLP.maxH+1,No_Countries);

for hh = 0:settingsLP.maxH
    
    yy_tplus_h = mat_single_leadlag(yy,hh) - mat_single_leadlag(yy,-1);  
    vec_yy_tplus_h = yy_tplus_h(:);
    
    vec_f_it = Countries_f_it;

    vec_Instr =[]; 
    for jjii = 1:size(Countries_Instr_f_it,3)
    vec_Instr = [vec_Instr squeeze(Countries_Instr_f_it(:,:,jjii))];
    end

    vec_Controls_t = [];
    for lag = 1:settingsLP.MaxLPLagsOwn
        mat_elem_dy = DemeanVariables(mat_single_leadlag(dy,-lag)); 
        vec_Controls_t = [vec_Controls_t mat_elem_dy(:)];
    end
    
    for lag = 1:settingsLP.MaxLPLagsOther
        mat_elem_Control = mat_single_leadlag(dx,-lag); 
        vec_Controls_t = [vec_Controls_t mat_elem_Control(:) ];
    end
    
    
    Select = isfinite(sum([vec_yy_tplus_h, vec_f_it, vec_Controls_t, vec_Instr],2)); 
    
    % B) Select Country on which you can estimate the interaction effect
    % vec_Instr has No_Countries*Q columns (Q blocks of No_Countries each).
    % CountriesToUse: country i is kept if it has >= MinEvents across all Q blocks.
    MinEvents = 3;
    Q = size(vec_Instr,2) / No_Countries;  % number of instruments
    CountriesToUse = false(1, No_Countries);
    for q = 1:Q
        idx = (q-1)*No_Countries + (1:No_Countries);
        CountriesToUse = CountriesToUse | (sum(vec_Instr(Select, idx) ~= 0, 1) >= MinEvents);
    end
    AllCountries = [1:No_Countries]; 
    SelectCountries = AllCountries(CountriesToUse); 



    % C) ESTIMATE
    % Build ZZ by selecting the CountriesToUse columns from each of the Q blocks
    ZZ_cols = [];
    for q = 1:Q
        idx = (q-1)*No_Countries + (1:No_Countries);
        ZZ_cols = [ZZ_cols idx(CountriesToUse)];
    end
    ZZ = vec_Instr(Select, ZZ_cols);
    ZZ = ZZ(:, (sum((ZZ~=0),1) > 1));
    ivfe = ivpanel(vec_id(Select,:), vec_year(Select,:), vec_yy_tplus_h(Select,:),...
                [vec_f_it(Select,CountriesToUse) vec_Controls_t(Select,:)], ...
                [ZZ], 'fe', 'endog', [1:size(vec_f_it(:,CountriesToUse),2)]);
    
    
    InteractCoeff(hh+1,SelectCountries) = ivfe.coef([1:size(vec_f_it(:,CountriesToUse),2)]); 
    InteractCoeffSE(hh+1,SelectCountries) = ivfe.stderr([1:size(vec_f_it(:,CountriesToUse),2)]);  
    
end

InteractTerm = zeros(settingsLP.maxH+1,No_Countries);
CollectSDTunits =zeros(settingsLP.maxH+1,1);
for rr = 1:size(InteractTerm,1); for cc = 1:size(InteractTerm,2)
        if isfinite(InteractCoeff(rr,cc))
            InteractTerm(rr,cc) = InteractCoeff(rr,cc)-nanmean(InteractCoeff(rr,:)'); 
        end
        CollectMeanunits(rr,1) = nanmean(InteractCoeff(rr,:)');
        CollectSDTunits(rr,1) = nanstd(InteractCoeff(rr,:)');
end; end

Collect_LP.InteractTerm = InteractTerm;
Collect_LP.InteractTermNONdemeaned = InteractCoeff;
Collect_LP.CollectMeanunits = CollectMeanunits;
Collect_LP.CollectSDTunits = CollectSDTunits;

InteractTermALL = zeros(T,No_Countries,settingsLP.maxH+1); 
for jj = 1:size(InteractTerm,1) 
    InteractTermALL(:,:,jj) = repmat(InteractTerm(jj,:),T,1); 
end


%%

maxH = settingsLP.maxH; 
MaxLPLagsOwn = settingsLP.MaxLPLagsOwn;
MaxLPLagsOther = settingsLP.MaxLPLagsOther; 

yy_mat = YY; 


%% Organize Inputs 

vec_id = a_id(:); 
vec_year = a_time(:); 

xx = XX_mat; 
dx = MatrixDiff(xx); 
dx = dx/median(nanstd(dx)); %% Normalization
    
InstrumentSelect = Instr_mat;

dy = yy_mat - mat_single_leadlag(yy_mat,-1); 

%% Organize Output 

WhenStartIRF = 0; 

%% ESTIMATE KBO WITH INTERACTION 

Collect_LP.ATE = ones(maxH+1,1); 
Collect_LP.ATEStdErr = 0*ones(maxH+1,1);

Collect_LP.Mediat = ones(maxH+1,1); 
Collect_LP.MediatStdErr = 0*ones(maxH+1,1);
Collect_LP.MediatSignificant = NaN*ones(maxH+1,1);

Collect_LP.CoeffMediatDOWN = ones(maxH+1,1); 
Collect_LP.StdErrMediatDOWN = ones(maxH+1,1);
Collect_LP.CoeffMediatUP = ones(maxH+1,1); 
Collect_LP.StdErrMediatUP = ones(maxH+1,1);     

ChooseSignificanceLevel = 10/100;


for hh = WhenStartIRF:maxH
            yy_tplus_h = mat_single_leadlag(yy_mat,hh) - mat_single_leadlag(yy_mat,-1);  
            vec_yy_tplus_h = yy_tplus_h(:);

            vec_Instr = [];
            for ivij = 1:size(InstrumentSelect,3);
                InstrumentSelect_i = InstrumentSelect(:,:,ivij);
                vec_Instr = [vec_Instr InstrumentSelect_i(:)];
                clear InstrumentSelect_i
            end
            
            vec_X_t = [dx(:)]; 
            vec_Controls_t = [];
            vec_Controls1Lag_t = [];

            for lag = 1:MaxLPLagsOwn
                mat_elem_dy = DemeanVariables(mat_single_leadlag(dy,-lag)); 
                vec_Controls_t = [vec_Controls_t mat_elem_dy(:)];
                if lag ==1 
                    vec_Controls1Lag_t = [vec_Controls1Lag_t mat_elem_dy(:)];
                end
            end
                        

            for pick = 1:size(Controls,3)
                pickControl = DemeanVariables(Controls(:,:,pick));
                for lag = 1:MaxLPLagsOther
                    mat_elem_Control = mat_single_leadlag(pickControl,-lag); 
                    vec_Controls_t = [vec_Controls_t mat_elem_Control(:) ];

                    if lag ==1 
                    vec_Controls1Lag_t = [vec_Controls1Lag_t mat_elem_Control(:)];
                    end
                end
            end
    
            if sum(size(DetermControls))~=0
            for pick = 1:size(DetermControls,3)
                if sum(sum(DetermControls(1:end-hh,:,pick)))>0 %% THIS DROP DET CONTROLS THAT ARE REDUNDANT
                pickDetControl = DetermControls(:,:,pick);
                % pickDetControl = DemeanVariables(DetermControls(:,:,pick));
                vec_Controls_t = [vec_Controls_t pickDetControl(:)];
                else 
                    fprintf(['Dropped Deterministic Control <<' num2str(pick) '>>, at Horizon [[' num2str(hh) ']]\n']); 
                end
            end
            end

            
                vec_ControlInteract_t = vec_Controls1Lag_t.*repmat(vec_X_t,1,size(vec_Controls1Lag_t,2));
                vec_InstrInteract_t = [];
                for instr_i = 1:size(vec_Instr,2);
                vec_InstrInteract_t = [vec_InstrInteract_t vec_Controls1Lag_t.*repmat(vec_Instr(:,instr_i),1,size(vec_Controls1Lag_t,2))];
                end
                

            switch KBOType.CaseInteraction
                case 'No_Interactions' 
                        Select = isfinite(sum([vec_yy_tplus_h, vec_X_t, vec_Controls_t, vec_Instr],2)); 
                    ivfe = ivpanel(vec_id(Select,:), vec_year(Select,:), vec_yy_tplus_h(Select,:),...
                        [vec_X_t(Select,:) vec_Controls_t(Select,:)], ...
                        [vec_Instr(Select,:)], 'fe', 'endog', 1);
                case 'With_Interactions'
                    Mediator_h = squeeze(InteractTermALL(:,:,hh+1));
                    vec_MediateInteract_t = vec_X_t.*Mediator_h(:);
                    vec_MediateInstrInteract_t = vec_Instr(:,1).*Mediator_h(:);
                    for ijij = 2:size(vec_Instr,2)
                    vec_MediateInstrInteract_t = [vec_MediateInstrInteract_t vec_Instr(:,ijij).*Mediator_h(:)];
                    end

                    Select = isfinite(sum([vec_yy_tplus_h, vec_X_t, vec_Controls_t, vec_ControlInteract_t, vec_Instr, vec_InstrInteract_t vec_MediateInstrInteract_t],2)); 
                    XX_t = [vec_X_t(Select,:) vec_MediateInteract_t(Select,:) vec_ControlInteract_t(Select,:) vec_Controls_t(Select,:)]; 
                    ZZ_t = [vec_Instr(Select,:) vec_MediateInstrInteract_t(Select,:) vec_InstrInteract_t(Select,:)];
                    ivfe = ivpanel(vec_id(Select,:), vec_year(Select,:), vec_yy_tplus_h(Select,:),...
                        XX_t, ZZ_t, 'fe', 'endog', [1:1+1+size(vec_Controls1Lag_t,2)]);            

            end


            
            Collect_LP.ATE(hh+1,1) = ivfe.coef(1,1); 
            Collect_LP.ATEStdErr(hh+1,1) = ivfe.stderr(1,1); 

            Collect_LP.Mediat(hh+1,1) = ivfe.coef(2,1); 
            Collect_LP.MediatStdErr(hh+1,1) = ivfe.stderr(2,1); 
            

            %% INTERACTED RESPONSES 
            
            StdSteps = 1*CollectSDTunits(hh+1); 
            Vector = [1 -StdSteps];                
            Collect_LP.CoeffMediatDOWN(hh+1,1) = Vector*ivfe.coef(1:2,1);
            Collect_LP.StdErrMediatDOWN(hh+1,1) = sqrt(Vector*ivfe.varcoef(1:2,1:2)*Vector');

            Vector = [1 StdSteps];                
            Collect_LP.CoeffMediatUP(hh+1,1) = Vector*ivfe.coef(1:2,1);
            Collect_LP.StdErrMediatUP(hh+1,1) = sqrt(Vector*ivfe.varcoef(1:2,1:2)*Vector');
            
            if hh == 0
                Collect_LP.delta_hat = ivfe.first_stage_coef(end, 1);
            end
            clear ivfe


        
end

FindSignificant = (abs(Collect_LP.Mediat(:,1)./Collect_LP.MediatStdErr(:,1))>=norminv(1-(ChooseSignificanceLevel/2))); 
Collect_LP.MediatSignificant(FindSignificant,1) = 1; 
