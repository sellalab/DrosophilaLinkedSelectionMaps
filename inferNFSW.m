
function [calc, preLHstats]  = inferNFSW( GEs, fdata, cfg_inf, bounds, masks, preLHstats, invar_files)


% cfg_inf.fmincon_retries = 0;

calc.init_params = bounds.init; % use the values just loaded in the previous script as the parameters here
calc.rangeL = bounds.rangeL; % use the L and H bounds just loaded in the previous script
calc.rangeH = bounds.rangeH;
calc.rangeL(bounds.fixed==1) = bounds.init(bounds.fixed==1); % copy the fixed params from L and H into the same indices
calc.rangeH(bounds.fixed==1) = bounds.init(bounds.fixed==1);

calc.DFs = sum( calc.rangeL<calc.rangeH ); % degrees of freedom i.e. total non-fixed params


iC = cfg_inf.chromosomes;

%% set algorithms

% just a formality of naming the algorithms and setting contraints based on
% what is set in default configuration


if isfield(cfg_inf,'constraint_u') && (cfg_inf.constraint_u == 1) % these are different contraints on u_del which can be set 
  f_constraint = @constraint_u_sum;
elseif isfield(cfg_inf,'constraint_u') && (cfg_inf.constraint_u == 2)
  f_constraint = @constraint_u_sum_ineq;
else
  f_constraint = []; % this is the setting where u_del has no contraints
end

if ~isfield(cfg_inf, 'opt_type') | isempty(cfg_inf.opt_type)
  cfg_inf.opt_type = {'active-set', 'interior-point', 'sqp'};
end

for i=1:length(cfg_inf.opt_type)
  opts{i} = optimset('Algorithm',cfg_inf.opt_type{i},     'Display','iter');
end

%% precalculate stats

if isempty(preLHstats)
  for c=iC
%     if isfield(cfg_inf, 'predef_idx_train') & ~isempty(cfg_inf.predef_idx_train)
    if exist('masks') & ~isempty(masks) % if masked are used, they are applied here (not used with the human data)
      predef_idx = masks{c};  
    else
      predef_idx = [];  %  now i have turned off masks entirely, skipping this step and (hopefully) saving some loading time
    end
    
    if ~isfield(GEs, 'BSbase')
      cur_BSbase = {[]};
    else
      cur_BSbase = GEs.BSbase(c,:); % struct with ranges over which BS coef is taking effect (Lj), size of effect (Bj) and config
    end
    
    if ~isfield(GEs, 'SWbase')
      cur_SWbase = {[]};
    else
      cur_SWbase = GEs.SWbase(c,:); % this contains coalescent rates due to sweeps on different s valued substitutions
      cur_gFocGrid = GEs.gFocGrid{c}; % this contains the postitions corresponding to the coalescent rates in the last struct
    end
    
    % edit for MEMORY MAP capability:
    % add "invar_files" as an input to use memory maps inside function
%      preLHstats{c} = calcPreLHStatsSw_memMap( c, cur_SWbase, cur_gFocGrid, cur_BSbase, fdata{c}, cfg_inf, invar_files  ); % predef_idx DM edit: remove weights (getting rid of masks)
%      sprintf('done with chr%s', c)
      preLHstats{c}    = calcPreLHStatsSw( cur_SWbase, cur_gFocGrid, cur_BSbase, fdata{c}, cfg_inf); % predef_idx DM edit: remove weights (getting rid of masks)
   
    
    % preLHstats contains gSWj: this is a matrix of A x N dimensions where
    % A is the number of selection coefficients (s) used and N is the number
    % of neutral sites in the dataset. Each N-length row has the predicted
    % coalscent rate from the .sw precalculated files
    % gBSj: this is a similarly A x N matrix where A is the the number of
    % background selection coefficients (t) and N is similarly the size of
    % the dataset. Each N-length row in this case has the predicted
    % decrease in diversity due to BS as a decimal i.e. -1% = -0.01

    %     calc.EgMutDiv(c) = preLHstats{c}.EgMutDiv;
    
end

%% initial values for the likelihood:

calc.init_LH = logCL_SW(preLHstats, calc.init_params, cfg_inf);
% calc.init_LH = logCL_SW_memMap(fdata, invar_files, preLHstats, calc.init_params, cfg_inf);

calc.best_iter = -1;
LH_best = calc.init_LH;
% params_temp  = calc.init_params;
% LH_temp      = calc.init_LH;


%% perform optimization:

opts_default = optimoptions('fmincon');

for k=1:length(opts)
  
  if calc.best_iter==-1
    calc.params(k,:)  = calc.init_params;
  else
    calc.params(k,:)  = calc.params(calc.best_iter,:);
  end
  
  retries = cfg_inf.fmincon_retries;
  calc.exitflag(k)    = 0;
  opts{k}.MaxIter     = opts_default.MaxIter;
  opts{k}.MaxFunEvals = 100*size(calc.params,2); %opts_default.MaxFunEvals;
  
  while ( retries >= 0 ) & ( calc.exitflag(k) == 0 )
    
    full_params = calc.params(k,:);
    ivariables  = find(bounds.fixed==0);
  
    if strcmp(cfg_inf.opt_type{k}, 'interior-point')
      full_params(ivariables) = epsilonFromBoundary( full_params(ivariables), calc.rangeL(ivariables), calc.rangeH(ivariables), 10^-9 );
    end
 
          % MODIFIED FOR MEMMAP: logCL_SW_memMap

    if isfield(cfg_inf,'constraint_u') && (cfg_inf.constraint_u == 1)
      [tvariables, calc.nLogLH(k), calc.exitflag(k), calc.opt{k}] = fmincon(@(variables)logCL_SW(preLHstats, full_params, cfg_inf, variables, ivariables), full_params(ivariables), [],[],[],[], calc.rangeL(ivariables), calc.rangeH(ivariables), @(variables)constraint_u_sum(full_params,variables,ivariables), opts{k});
%       [tvariables, calc.nLogLH(k), calc.exitflag(k), calc.opt{k}] = fmincon(@(variables)logCL_SW_memMap(fdata, invar_files, preLHstats, full_params, cfg_inf, variables, ivariables), full_params(ivariables), [],[],[],[], calc.rangeL(ivariables), calc.rangeH(ivariables), @(variables)constraint_u_sum(full_params,variables,ivariables), opts{k});
    elseif isfield(cfg_inf,'constraint_u') && (cfg_inf.constraint_u == 2)
      [tvariables, calc.nLogLH(k), calc.exitflag(k), calc.opt{k}] = fmincon(@(variables)logCL_SW(preLHstats, full_params, cfg_inf, variables, ivariables), full_params(ivariables), [],[],[],[], calc.rangeL(ivariables), calc.rangeH(ivariables), @(variables)constraint_u_sum_ineq(full_params,variables,ivariables), opts{k});
%       [tvariables, calc.nLogLH(k), calc.exitflag(k), calc.opt{k}] = fmincon(@(variables)logCL_SW_memMap(fdata, invar_files, preLHstats, full_params, cfg_inf, variables, ivariables), full_params(ivariables), [],[],[],[], calc.rangeL(ivariables), calc.rangeH(ivariables), @(variables)constraint_u_sum_ineq(full_params,variables,ivariables), opts{k});
    
    else
        
      % THE CASE WE CURRENTLY USE: NO RESTRICTION PUT ON U DELETERIOUS...
      [tvariables, calc.nLogLH(k), calc.exitflag(k), calc.opt{k}] = fmincon(@(variables)logCL_SW(preLHstats, full_params, cfg_inf, variables, ivariables), full_params(ivariables), [],[],[],[], calc.rangeL(ivariables), calc.rangeH(ivariables), [], opts{k});
%       [tvariables, calc.nLogLH(k), calc.exitflag(k), calc.opt{k}] = fmincon(@(variables)logCL_SW_memMap(fdata, invar_files, preLHstats, full_params, cfg_inf, variables, ivariables), full_params(ivariables), [],[],[],[], calc.rangeL(ivariables), calc.rangeH(ivariables), [], opts{k});

    end
    calc.params(k,ivariables) = tvariables;
    
    if calc.exitflag(k) == 0
      opts{k}.MaxFunEvals = 2*opts{k}.MaxFunEvals; %[num2str(2*str2num(opts{k}.MaxFunEvals(1:find(opts{k}.MaxFunEvals=='*')-1))) '*numberOfVariables'];
      opts{k}.MaxIter     = 2*opts{k}.MaxIter;
    end
    retries = retries-1;
  end
  
  if ( LH_best >= calc.nLogLH(k) ) & ( calc.exitflag(k) >= 0 )
    calc.best_iter = k;
    LH_best = calc.nLogLH(k);
  end
  %   if LH_temp >= calc.nLogLH(k)  &  calc.exitflag(k) > 0
  % %   if LH_temp >= calc.nLogLH(k)
  %     params_temp = calc.params{k};
  %     LH_temp     = calc.nLogLH(k);
  %     calc.best_iter = k;
  %   end
  
end

%%

calc.nLogLH = calc.nLogLH / 10^5; % remove the 10^5 factor from inside logCL_SW

if calc.best_iter~=-1
  calc.best_iter = find( calc.nLogLH == min(calc.nLogLH), 1, 'last' );
  [~, calc.samples, calc.LHFuncStats] = logCL_SW(preLHstats, calc.params(calc.best_iter,:), cfg_inf);
%   [~, calc.samples, calc.LHFuncStats] = logCL_SW_memMap(fdata, invar_files, preLHstats, calc.params(calc.best_iter,:), cfg_inf);
  calc.logLH = -calc.nLogLH( calc.best_iter );
else
  calc.samples = [0 0 0];
  calc.logLH = 0;
end



% if isfield(cfg_inf, 'predef_idx_train')
%   % this is a patch in order to save space - theoreticaly we need to track
%   % the exact list of codons used for the inference, but we settle with
%   % keeping this field only at the 'anals_obs' struct for all corresponding
%   % inference runs
%   calc.config = rmfield(cfg_inf, 'predef_idx_train');
% else
%   bigo = 7;
% end
% % if isfield(calc.config, 'predef_idx_test')
% %   calc.config = rmfield(calc.config, 'predef_idx_test');
% % end
% % if isfield(calc.config, 'pos')
% %   calc.config = rmfield(calc.config, 'pos');
% % end



calc.config.inf = cfg_inf;
calc.config.GEs = GEs.cfg;


end

