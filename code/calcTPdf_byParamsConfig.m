
function [v_t, w_t, w_t_rel, u_del_max, u_del, u_del_rel, nsites_del] = calcTPdf_byParamsConfig( t_params, BSbase_params, config )

global MLParamsStruct;

v_t = [];
w_t = [];
w_t_rel = [];
u_del = [];

minimal_log10_t = -10 + 0.001; % reset this to be just above the null -10 that is used on fixed params...

u_del_max = config.u_del_max; % just take this from the original value

nsites_del = [];
% BSbase_params = number of annos (i.e. exonic, nonexonic, etc)
for a=1:length(BSbase_params)
  % number of sites del is just all sites in a given anno since they are all subject to possible mutations
  nsites_del(a) = BSbase_params{a}.nsites; 
end

for a=1:length(BSbase_params)

  % BSanno2param_mapping is simply the current anno mapping ([1 2] for human)
  aa = config.BSanno2param_mapping(a); 
  % v_t{a} = the vector of t vals
  v_t{a} = BSbase_params{a}.t; 
  nvectors = length(v_t{a});
  % confirm that there are not too many precalc grid points for some reason
  assert( nvectors <= MLParamsStruct.bsparam_masses ,'precalculated base contains more vectors than possible degrees of freedom' ); 
  % for each set of t params get the current weights based on LHM. this only considers the weights that we have not fixed
  cur_weights = t_params((aa-1)*MLParamsStruct.bsparam_annolen+[1:nvectors]); 
  irangeL = find(cur_weights <= minimal_log10_t); % find weights that have gone below the minimum (if any)
  absw_t = 10.^cur_weights; % convert the weight out of log to a base 10 value
  absw_t(irangeL) = 0; % points that have dropped below the cutoff are set to 0 likelihood
  
  if sum(absw_t)>0 % this evaluates true until all of the weights have dropped below the cutoff minimal_log10_t
    w_t{a} = absw_t/sum(absw_t); % weights are distributed over those remaining above the threshold
  else
    w_t{a} = zeros(size(absw_t));
  end
  
  u_del(a)     = sum(absw_t); % deleterious mut rate set at sum of the weights (?)
  u_del_rel(a) = u_del(a) / u_del_max; % relative mutation rate to the default max (should we change this in humans?)
  
  w_t_rel{a} = absw_t ./ BSbase_params{a}.u; % relative WEIGHTS are given by ratio of weight to bsparam u
end

% end

end

