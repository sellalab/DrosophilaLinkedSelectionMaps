
function [neg_log_P, neg_log_samples, stats, vlogL] = logCL_SW(preCalc, full_params, config, variables, ivariables)

global MLParamsStruct;

% To prevent zero probability of observing a SNP (which would ruin the
% maximization) we put a lower bound on the relative predicted diversity to
% be 0.1% of the maximal value (=1).
minRed = 0.001; 

calc_stats = (nargout > 2);


if nargin < 5
  variables  = [];
  ivariables = [];
end

params             = full_params;
params(ivariables) = variables;


iC = config.chromosomes;


TAU_div     = params(3); % pi0

% A patch to push the maximization algorithm out of regions in the
% parameter space where it estimates for some reason some of the parameters
% as NaN.
if sum(isnan(params))
  neg_log_P = 100000000000000000000;
  return;
end

neg_log_P = 0;
neg_log_samples = [0 0 0];

%% predict diversity reduction based on all params still in play

for c=iC
  % (diversity reduction prediction):
  % the first argument refers to the flag only_calc_Red, which, when turned
  % on, skips two individual calculations of SW and BS. Keep on for now:

  DivRedPred = composeLSMapFromElements( 1, preCalc{c}.gSWj, preCalc{c}.gBSj, params, preCalc{c}.SWparams, preCalc{c}.BSparams, preCalc{c}.EgMutDiv, config ); % DM edit: set only_calc_Red to 0 for BS or SW only maps
  
  theta0   = preCalc{c}.gMutDiv  / TAU_div; % proxy for 4*Ne*u here?
  
  pi0 = max( theta0./(1+theta0), 10^-5); % neutral pi = theta / 1+ theta (theta should only change if local u changes, currently fixed for human)
  
  pii        = min( pi0 .* max( DivRedPred.Red, minRed ), 0.99 ); % update the current observed pi based on the diversity reduction factor that has been calculated
  
  vlogL{c}      = double(preCalc{c}.samplesHom).*log(1-pii) + double(preCalc{c}.samplesHet).*log(pii);
  cneg_log_P{c}  = -sum(vlogL{c});
  
  cneg_log_samples{c} = [preCalc{c}.samplesHomSum+preCalc{c}.samplesHetSum  preCalc{c}.samplesHomSum  preCalc{c}.samplesHetSum];  % [all_sites homo_sites het_sites]
  csites{c} = length(preCalc{c}.samplesHet);
  
  if nargout>2 % only taking 1 arg out in infer_NFSW
    SumMutProx(c) = preCalc{c}.SumMutProx;
    nsamples(c) = preCalc{c}.samplesSum;
    nhetpairs(c)= preCalc{c}.samplesHetSum;
  end
end


for c=iC
  neg_log_samples = neg_log_samples + cneg_log_samples{c};
  neg_log_P       = neg_log_P       + cneg_log_P{c};
end

% A patch to prevent too low likelihood values, for numeric reasons. This
% patch is complemented by a corresponding division by 10^5 at the
% encapsulating level.
neg_log_P = neg_log_P/neg_log_samples(1) * 10^5;

% if isnan(neg_log_P) | isinf(neg_log_P)
%   bigo = 7;
% end


% [TAU_div params(59:61) neg_log_P]  % for minimal BS params
% [TAU_div params(11:2:15) neg_log_P]  % for minimal SW params

% [params(1) params(3) params(7:end) neg_log_P/1000000] %;params(1) params(3) w_s neg_log_P/1000000
% [params(1:3) neg_log_P/1000000]

%%

if nargout>2 % only taking one argument out in infer_NFSW
  stats.sites       = csites;
  stats.samplepairs = cneg_log_samples;
  stats.nlogLH      = cneg_log_P;
  stats.EMutProx    = sum(SumMutProx)  / sum(nsamples);
  stats.EHet        = sum(nhetpairs) / sum(nsamples);
end


bigo = 7;

end
