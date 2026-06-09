function [Collect_LP]=Panel_LP(YY,XX,Controls,DetermControls,a_id,a_time,settingsLP, interaction)
% Inputs:   YY: TxN dependent variable (panel outcome)
%           XX: TxN matrix of regressors OR Tx1 when common across units
%               (interpreted inside as shock variable whose first-difference is used)
%           Controls: TxNxK 3D matrix of additional controls (panel)
%           DetermControls: TxNxKd 3D matrix of deterministic controls (panel),
%                          can be [] if none
%           a_id:   TxN matrix of unit identifiers (panel structure)
%           a_time: TxN matrix of time identifiers (panel structure)
%           settingsLP: struct with LP settings:
%               settingsLP.maxH          = maximum horizon H (integer)
%               settingsLP.MaxLPLagsOwn  = number of lags of ΔY used as controls (integer)
%               settingsLP.MaxLPLagsOther= number of lags of other controls (integer)
%           intercation: 'No_Interactions' or 'With_Interactions'
%
% Outputs:  Collect_LP = struct()
%           Collect_LP.Coeff   = (H+1)x1 IRF coefficients (baseline shock coefficient)
%           Collect_LP.StdErr  = (H+1)x1 standard errors of baseline IRF coefficients
%
%           Collect_LP.CoeffInteractionUP   = (H+1)xNoInteractions IRFs evaluated at
%                                            +Magnitude for each interaction term
%           Collect_LP.StdErrInteractionUP  = (H+1)xNoInteractions standard errors for UP IRFs
%           Collect_LP.CoeffInteractionDOWN = (H+1)xNoInteractions IRFs evaluated at
%                                            -Magnitude for each interaction term
%           Collect_LP.StdErrInteractionDOWN= (H+1)xNoInteractions standard errors for DOWN IRFs
%

%% SPECIFY TYPE 

KBOType.CaseInteraction = interaction; 

%%

NoInteraction = size(Controls,3)+1;
maxH = settingsLP.maxH; 
MaxLPLagsOwn = settingsLP.MaxLPLagsOwn;
MaxLPLagsOther = settingsLP.MaxLPLagsOther; 

yy_mat = YY; 
[~,NN] = size(yy_mat);

%% If XX Common Across Countries -- Create Panel 
if size(XX,2) == 1
    XX_mat = repmat(XX,1,NN); 
else 
    XX_mat = XX; 
end


%% Organize Inputs 

vec_id = a_id(:); 
vec_year = a_time(:); 

xx = XX_mat; 
dx = xx; 
dx = dx/median(nanstd(dx)); %% Normalization

dy = yy_mat - mat_single_leadlag(yy_mat,-1); 

%% Organize Output 

WhenStartIRF = 0; 
Collect_LP.Coeff = ones(maxH+1,1); 
Collect_LP.StdErr = 0*ones(maxH+1,1);

Collect_LP.CoeffInteractionUP = ones(maxH+1,NoInteraction); 
Collect_LP.StdErrInteractionUP = 0*ones(maxH+1,NoInteraction);
Collect_LP.CoeffInteractionDOWN = ones(maxH+1,NoInteraction); 
Collect_LP.StdErrInteractionDOWN = 0*ones(maxH+1,NoInteraction);

for hh = WhenStartIRF:maxH
            yy_tplus_h = mat_single_leadlag(yy_mat,hh) - mat_single_leadlag(yy_mat,-1);  
            vec_yy_tplus_h = yy_tplus_h(:);
            
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

            switch KBOType.CaseInteraction
                case 'No_Interactions' 
                        Select = isfinite(sum([vec_yy_tplus_h, vec_X_t, vec_Controls_t],2)); 
                    olsfe = panel_ols_fe(vec_id(Select,:), vec_year(Select,:), vec_yy_tplus_h(Select,:),...
                        [vec_X_t(Select,:) vec_Controls_t(Select,:)]);
                case 'With_Interactions' 
                    Select = isfinite(sum([vec_yy_tplus_h, vec_X_t, vec_Controls_t, vec_ControlInteract_t],2)); 
                    olsfe = panel_ols_fe(vec_id(Select,:), vec_year(Select,:), vec_yy_tplus_h(Select,:),...
                        [vec_X_t(Select,:) vec_ControlInteract_t(Select,:) vec_Controls_t(Select,:)]);            
            end


            Collect_LP.Coeff(hh+1,1) = olsfe.coef(1,1); 
            Collect_LP.StdErr(hh+1,1) = olsfe.stderr(1,1); 
            
            
            switch KBOType.CaseInteraction
                case 'With_Interactions' 
                    for Interact = 1:NoInteraction
                        if Interact == 1
                            MagnitInteract = median(std(dy,'omitnan'));
                        else
                            MagnitInteract = median(std(Controls(:,:,Interact-1),'omitnan'));
                        end
                        
                        Vector = zeros(1,NoInteraction); 
                        Vector(Interact) = MagnitInteract;
                        Vector = [1 Vector];
                        
                        Collect_LP.CoeffInteractionUP(hh+1,Interact) = Vector*olsfe.coef(1:1+NoInteraction,1);
                        Collect_LP.StdErrInteractionUP(hh+1,Interact) = sqrt(Vector*olsfe.varcoef(1:1+NoInteraction,1:1+NoInteraction)*Vector');
                        
                        Vector = zeros(1,NoInteraction); 
                        Vector(Interact) = -MagnitInteract;
                        Vector = [1 Vector];
                        
                        Collect_LP.CoeffInteractionDOWN(hh+1,Interact) = Vector*olsfe.coef(1:1+NoInteraction,1);
                        Collect_LP.StdErrInteractionDOWN(hh+1,Interact) = sqrt(Vector*olsfe.varcoef(1:1+NoInteraction,1:1+NoInteraction)*Vector');
                    
                    end
                end
            clear olsfe
        
end        

end

function [olsfe]=panel_ols_fe(id, time, y, X)

[isBalanced, idx, n] = isbalancedpanel(id, time);

id = id(idx);
y = y(idx,:);
X = X(idx,:);

N = size(y,1);
k = size(X,2);

ytr = y - groupmeans(id,y,'replicate',1);
Xtr = X - groupmeans(id,X,'replicate',1);

coef = (Xtr'*Xtr)\Xtr'*ytr;
res = ytr - Xtr*coef;
resdf = (N-n)-k;
resvar = (res'*res) ./ resdf;
varcoef = resvar * ((Xtr'*Xtr)\eye(k));
stderr = sqrt(diag(varcoef));

olsfe.coef = coef;
olsfe.stderr = stderr;
olsfe.varcoef = varcoef;

end
