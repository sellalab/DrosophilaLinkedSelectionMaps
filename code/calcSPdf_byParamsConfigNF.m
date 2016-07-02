
function [v_s, coalrate_s] = calcSPdf_byParamsConfigNF(s_params, SWbase_params, config)

global MLParamsStruct;

v_s = [];
coalrate_s = [];


for a=1:length(SWbase_params)
  v_s{a} = SWbase_params{a}.s; % v_s is just the set of s params
  coalrate_s{a} = zeros(size(v_s{a})); % coalrate initialized as a vector of zeros the size of the vector of s

  nvectors = length(v_s{a});
  assert( nvectors <= MLParamsStruct.swparam_masses ,'precalculated base contains more vectors than possible degrees of freedom' );

  % for each annotation, find the appropriate parameters set accodring to the mapping from the configutaion, iSWanno
  if ~isempty(config.SWanno2param_mapping)
    aa = config.SWanno2param_mapping(a); % aa just equals the index of the parameter mapping (i.e. 1 and 2 in the simple case we're starting with)
    coalrate_s{a} = max(s_params((aa-1)*(MLParamsStruct.swparam_masses+1)+[1:nvectors]),0); 
  else 
    coalrate_s{a} = 0;
  end

end


end
