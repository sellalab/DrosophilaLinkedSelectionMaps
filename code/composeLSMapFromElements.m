
function DivRedPred = composeLSMapFromElements( only_calc_Red, gSWj, gBSj, params, SWbase_params, BSbase_params, EgMutDiv, config )

global MLParamsStruct;

if isempty(gSWj) & isempty(gBSj)
  DivRedPred = [];
  return;
end

if ~isempty(gSWj)
  L = size(gSWj{1},2);
else
  L = size(gBSj{1},2);
end


% integrate params from both GEs and inference results
DivRedPred.params = integrateParams_Inference_N_GEs( params, config, BSbase_params, SWbase_params, 0);

% prepare mutation rate variation proxy
 eTheta0 = EgMutDiv / DivRedPred.params.tau_div;

%% create basic maps for sweeps and bs

% compose the basic map of coalescent due to SW - the minimal calculations possible
cSW  = zeros([L 1]); %
for a=1:length(gSWj)
  if sum( DivRedPred.params.SW.coalrate_s{a}(1:size(gSWj{a},1)) ) > 0
    cSW         = cSW + gSWj{a}' * DivRedPred.params.SW.coalrate_s{a}';%(2:2:10)';
  end
end

% compose the basic map of coalescent due to BS - the minimal calculations possible
cBS  = zeros([L 1]);
% length(gBSj) is the # of pre-calc maps used, also same as the number of selection coefs
for a=1:length(gBSj)
  % size(gBSj{a},1) is the number of sites, which mirrors number of neutral sites, per pre-calc map
  if sum( DivRedPred.params.BS.w_t_rel{a}(1:size(gBSj{a},1)) ) > 0
    % tranpose so that inner dimensions agree (i.e. # weights should equal # selection coefs used)
    cBS = cBS + gBSj{a}' * DivRedPred.params.BS.w_t_rel{a}';
  end
end


%% create individual maps for each effect in isolation (SW and BS)
% USED IN FINAL RESULT PRINT OUT
% FOR CL CALCULATION, EFFECTS ARE COMBINED
% compose additional maps of coalescent due to SW, by annotation and fitness effect, and also component-specific diversity maps

DivRedPred.cSW = cSW;

if ~only_calc_Red & isfield(DivRedPred.params, 'SW') % if only_calc_Red flag is NOT set (i.e. 0), calculate the individual components (?)
  
  DivRedPred.cSWs = zeros([L length(DivRedPred.params.SW.w_s{1})]);
  
  for a=1:length(gSWj)
    
    if sum( DivRedPred.params.SW.coalrate_s{a}(1:size(gSWj{a},1)) ) > 0
      cSWan{a}    = gSWj{a}' * DivRedPred.params.SW.coalrate_s{a}';
    else
      cSWan{a}    = zeros([L 1]);
    end
    
    for k=1:size(gSWj{a},1)
      DivRedPred.SWj{a}(:,k) = (1 + eTheta0) ./ (1 +  gSWj{a}(k,:)'*DivRedPred.params.SW.coalrate_s{a}(k) + eTheta0 );
      DivRedPred.cSWs(:,k)   = DivRedPred.cSWs(:,k) + gSWj{a}(k,:)'*DivRedPred.params.SW.coalrate_s{a}(k);
    end
    DivRedPred.cSWan{a}    = cSWan{a};
    DivRedPred.SWan{a}     = (1 + eTheta0) ./ (1 + cSWan{a} + eTheta0);
    
    % calculate derivatives of the parameters - the average impact of LS on
    % diversity across neutral sites
    for k=1:size(gSWj{a},1)
      DivRedPred.params.SW.eRelRed{a}(k)   = 1 - nanmean( DivRedPred.SWj{a}(:,k) );
      DivRedPred.params.SW.eRelCoalrate{a}(k)   = nanmean( gSWj{a}(k,:) ) * DivRedPred.params.SW.coalrate_s{a}(k);
    end
    DivRedPred.params.SW.eRelRedAnno(a)      = 1 - nanmean( DivRedPred.SWan{a} );
    DivRedPred.params.SW.eRelCoalrateAnno(a) = sum( DivRedPred.params.SW.eRelCoalrate{a}(k) );
    
  end
  
  DivRedPred.SW         = (1 + eTheta0) ./ (1 + DivRedPred.cSW  + eTheta0);
  DivRedPred.SWs        = (1 + eTheta0) ./ (1 + DivRedPred.cSWs + eTheta0);
  
end


% compose additional maps of coalescent due to BS, by annotation and fitness effect, and also component-specific diversity maps

DivRedPred.BS         = exp(cBS); % exp(effect of bs)

if ~only_calc_Red & isfield(DivRedPred.params, 'BS')
  DivRedPred.BSt = ones([L length(DivRedPred.params.BS.w_t{1})]);
  
  for a=1:length(gBSj)
    
    if sum( DivRedPred.params.BS.w_t_rel{a}(1:size(gBSj{a},1)) ) > 0
      cBSan{a}    = gBSj{a}' * DivRedPred.params.BS.w_t_rel{a}'; % this mirrors the above calc for cBS
    else
      cBSan{a}    = zeros([L 1]);
    end
    
    for k=1:size(gBSj{a},1) %find(DivRedPred.params.BS.w_t{a}>0)
      if DivRedPred.params.BS.w_t{a}(k)>0
        DivRedPred.BSj{a}(:,k) = exp(DivRedPred.params.BS.w_t_rel{a}(k) * gBSj{a}(k,:)');
        %         DivRedPred.BSj{a}(:,k) = exp((DivRedPred.params.BS.u_del{a}/BSbase_params{a}.u(k)*DivRedPred.params.BS.w_t{a}(k)) * gBSj{a}(k,:)');
        DivRedPred.BSt(:,k)  = DivRedPred.BSt(:,k).* DivRedPred.BSj{a}(:,k);
      else
        DivRedPred.BSj{a}(:,k) = ones([L 1]);
      end
    end
    DivRedPred.BSan{a}   = exp(cBSan{a}); % exponentiate the result for these weights to incorporate into the prediction for diversity reduction
    
    % calculate derivatives of the parameters - the average impact of LS on
    % diversity across neutral sites    
    for k=1:size(gBSj{a},1) %find(DivRedPred.params.BS.w_t{a}>0)
      if DivRedPred.params.BS.w_t{a}(k)>0
        DivRedPred.params.BS.eRelRed{a}(k)   = 1 - nanmean( (1 + eTheta0) ./ (DivRedPred.BSj{a}(:,k).^-1 + eTheta0) ); % expected reduction in diversity after the effects of sweeps (hence 1 - ... )
        DivRedPred.params.BS.eRelCoalrate{a}(k)   = nanmean( DivRedPred.BSj{a}(:,k).^-1 ) - 1;
        %         DivRedPred.params.BS.eRelCoalrate{a}(k)   = 1 - nanmean( DivRedPred.BSj{a}(:,k) );
      else
        DivRedPred.params.BS.eRelRed{a}(k)   = 0;
        DivRedPred.params.BS.eRelCoalrate{a}(k)  = 0;
      end
    end
    DivRedPred.params.BS.eRelRedAnno(a) = 1 - nanmean( (1 + eTheta0) ./ (DivRedPred.BSan{a}.^-1 + eTheta0) );
    DivRedPred.params.BS.eRelCoalrateAnno(a) = sum( DivRedPred.params.BS.eRelCoalrate{a} );
    
  end

%   DivRedPred.BS         = exp(cBS);

end


DivRedPred.Red   = ( 1 + eTheta0 ) ./ ( 1./DivRedPred.BS + DivRedPred.cSW + eTheta0 ); % this will mirror the rate from BS alone in my current settings




end

