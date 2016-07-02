
function [c,ceq] = constraint_u_sum( full_params, varialbes, ivariables )


global MLParamsStruct;


if nargin < 3
  variables = [];
  ivariables = [];
end

full_params(ivariables) = variables;


for k=1:MLParamsStruct.bsparam_annotations
  if full_params(MLParamsStruct.bsparam_imaxu) > -10;
    ceq(k)     = log10(sum(10.^full_params(MLParamsStruct.bsparam_imasses(k)+[0:MLParamsStruct.bsparam_masses-1]))) - full_params(MLParamsStruct.bsparam_imaxu(k));
  else
    ceq(k)     = 0;
  end
  c(k) = 0;
end

end