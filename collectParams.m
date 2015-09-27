
function param_stats = collectParams( calc, GEs, config )


param_stats.fit.samples  = calc.samples;
param_stats.fit.het      = calc.LHFuncStats.EHet;
param_stats.fit.EMutProx = calc.LHFuncStats.EMutProx;
param_stats.fit.logLH    = calc.logLH;
param_stats.fit.AIC      = 2*calc.DFs - 2*param_stats.fit.logLH*param_stats.fit.samples(1);






SWbase_params = [];
if ~isempty(GEs.SWbase{1})
  for a=1:size(GEs.SWbase,2)
    SWbase_params{a}.s    = GEs.SWbase{1,a}.cfg.FE_grid;
    SWbase_params{a}.Ne0  = GEs.SWbase{1,a}.cfg.Ne0; %SWbase{a}.params{k}.Ne0;
    SWbase_params{a}.remote2closeDiv = config.ratio_tau_mutprox_tau_sweeps;
  end
else
  SWbase_params{1}.s = 1;
  SWbase_params{1}.Ne0 = 1;
  SWbase_params{1}.remote2closeDiv = 1;
end

BSbase_params = [];
if ~isempty(GEs.BSbase{1})
  for a=1:size(GEs.BSbase,2)
    BSbase_params{a}.u         = GEs.BSbase{1,a}.cfg.u_del;
    BSbase_params{a}.t         = GEs.BSbase{1,a}.cfg.FE_grid;
    BSbase_params{a}.nsites    = GEs.BSbase{1,a}.cfg.anno_len;
  end
else
  BSbase_params{1}.u         = 10^-10;
  BSbase_params{1}.t         = 1;
  BSbase_params{1}.nsites    = 1;
end


if calc.best_iter~=-1
  param_stats.params = integrateParams_Inference_N_GEs( calc.params(calc.best_iter,:), config, BSbase_params, SWbase_params, 0 );
else
  param_stats.params = [];
end
% param_stats.params.tau = calc.params(3);
param_stats.params.pi0_by_d = 1/param_stats.params.tau_div;


end

