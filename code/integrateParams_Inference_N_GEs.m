
function params = integrateParams_Inference_N_GEs( Inf_params, Inf_config, BSbase_params, SWbase_params, isX, iFE_grid_X )

global MLParamsStruct;


params.tau_div       = Inf_params(MLParamsStruct.tau_pos);

%% BS param adjustments

[ params.BS.v_t,...
  params.BS.w_t,...
  params.BS.w_t_rel,...
  params.BS.u_del_max,...
  params.BS.u_del,...
  params.BS.u_del_rel,...
  params.BS.nsites_del] = calcTPdf_byParamsConfig(   Inf_params(MLParamsStruct.bsparam_offset+[1:MLParamsStruct.bsparam_annolen*MLParamsStruct.bsparam_annotations]), BSbase_params, Inf_config );
  params.BS.u_del_max   = Inf_config.u_del_max; % seems redundant, already set in the above function...
for a=1:length(params.BS.v_t) % length is currently six (max would be 11)
  params.BS.U_del(a)    = params.BS.u_del(a)*params.BS.nsites_del(a); % this will change as u_del changes over the maximization
  params.BS.Et(a)       =       sum(params.BS.w_t{a}.*params.BS.v_t{a} ); % mean t
  params.BS.St(a)       = sqrt( sum(params.BS.w_t{a}.*params.BS.v_t{a}.^2) - params.BS.Et(a)^2 ); % variance t
end

%% SW param adjustments

[ params.SW.v_s,...
  params.SW.coalrate_s] = calcSPdf_byParamsConfigNF( Inf_params(MLParamsStruct.swparam_offset+[1:MLParamsStruct.swparam_annolen*MLParamsStruct.swparam_annotations]), SWbase_params, Inf_config );
if ~isempty(SWbase_params)
  params.SW.Ne0               = SWbase_params{1}.Ne0;
  params.SW.remote2closeDiv   = SWbase_params{1}.remote2closeDiv; % this is the only use of remote2closeDiv in the main inference (that is not commented atm)
  params.SW.tau_subs          = params.tau_div / params.SW.remote2closeDiv; % used next in assignment here...
end

for a=1:length(params.SW.v_s)
  params.SW.alpha_over_tauS(a) = sum( params.SW.coalrate_s{a} ); % (total increase to coalescent rate induced by the sweep params) 
  if params.SW.alpha_over_tauS(a) > 0
    params.SW.w_s{a} = params.SW.coalrate_s{a} / params.SW.alpha_over_tauS(a); % weights are set relative to the total effect of sweeps on coalescent rate (tau_S)
  else
    params.SW.w_s{a} = params.SW.coalrate_s{a};
  end
  
  params.SW.alpha(a)    = params.SW.alpha_over_tauS(a) * params.SW.tau_subs; % downstream use of remote2closeDiv: sum of the coalrates * tau_subs (alpha = fraction of beneficial mutations?) a single fraction for each class (coding, non, etc)
  params.SW.alpha_s{a}  = params.SW.alpha(a)*params.SW.w_s{a}; % weights multiplied by alpha gives alpha_s: fraction of mutations with param s that are beneficial.
  
  params.SW.Es0(a)      =       sum(params.SW.alpha(a)*params.SW.w_s{a}.*params.SW.v_s{a});
  params.SW.Ss0(a)      = sqrt( sum(params.SW.alpha(a)*params.SW.w_s{a}.*params.SW.v_s{a}.^2) - params.SW.Es0(a)^2 );
  params.SW.Es(a)       =       sum(params.SW.w_s{a}.*params.SW.v_s{a} );
  params.SW.Ss(a)       = sqrt( sum(params.SW.w_s{a}.*params.SW.v_s{a}.^2) - params.SW.Es(a)^2 );  
end


%% optional for sex chromosomes

if nargin < 4
  isX = 0;
end
if nargin < 5
  iFE_grid_X = [1:MLParamsStruct.bsparam_masses];
end

if isX % this is false for now (X chromosome; scale things by 4/3, etc.)
  
  for a=1:length(params.SW.alpha_over_tauS)
%     params.SW.w_s_X{a} = zeros(size(params.SW.w_s{a}));
    grid_s = params.SW.v_s{a}(iFE_grid_X(end:-1:1));
    ws     = params.SW.w_s{a}(iFE_grid_X(end:-1:1));
    [wsX, vsX_log10] = resampleDiscreteDistribution( log10( 4/3 * grid_s ), ws, log10( grid_s ),  log10( [grid_s(1)/10 grid_s(end)*10] ) );
    %   [wsX, vsX_log10] = resampleDiscreteDistribution(        4/3 * grid_s  , ws,        grid_s  ,         [grid_s(1)/10 grid_s(end)*10]   );
    params.SW.w_s_X{a} = wsX(end:-1:1);
    params.SW.v_s_X{a} = 10.^vsX_log10(end:-1:1);
  end
  
  for a=1:length(params.BS.u_del)
%     params.BS.w_t_X{a} = zeros(size(params.BS.w_t{a}));
    grid_t = params.BS.v_t{a}(iFE_grid_X(end:-1:1));
    wt     = params.BS.w_t{a}(iFE_grid_X(end:-1:1));
    [wtX, vtX_log10] = resampleDiscreteDistribution( log10( 4/3 * grid_t ), wt, log10( grid_t ),  log10( [grid_t(1)/10 grid_t(end)*10] ) );
    params.BS.w_t_X{a} = wtX(end:-1:1);
    params.BS.v_t_X{a} = 10.^vtX_log10(end:-1:1);
  end
  
end





end

