
function gl_initMLParamsStruct()

global MLParamsStruct;

MLParamsStruct.tau_pos = 3;
MLParamsStruct.swparam_offset = 10; % DM what is the offset mean?
MLParamsStruct.swparam_masses = 11; % max number of s values used in the inference
MLParamsStruct.swparam_annotations = 4; % max number of annotations used, currently 4
MLParamsStruct.swparam_annolen = MLParamsStruct.swparam_masses + 1; % 11 positions in the vector for max number of s values, + 1 for ?
MLParamsStruct.bsparam_offset = MLParamsStruct.swparam_offset + MLParamsStruct.swparam_annolen*MLParamsStruct.swparam_annotations; % start bsparam at position 58 in the vector after the end of swparams
MLParamsStruct.bsparam_masses = 11; % same max number of s values
MLParamsStruct.bsparam_annotations = 5; %changed for human DM-- restored to 4
MLParamsStruct.bsparam_annolen = MLParamsStruct.bsparam_masses + 1; % same length of 12, 11 s and 1 extra something...
MLParamsStruct.length = MLParamsStruct.bsparam_offset + (MLParamsStruct.bsparam_masses+1)*MLParamsStruct.bsparam_annotations; % 48 BS + 48 SW params + 10 additional in the beginning

MLParamsStruct.swparam_imasses = MLParamsStruct.swparam_offset  + 1 + [0:MLParamsStruct.swparam_annotations-1]*MLParamsStruct.swparam_annolen; % indices of the start point of each subset of SW params per anno
MLParamsStruct.swparam_imaxu   = MLParamsStruct.swparam_imasses + MLParamsStruct.swparam_masses; % end points of each subset of SW params
MLParamsStruct.bsparam_imasses = MLParamsStruct.bsparam_offset  + 1 + [0:MLParamsStruct.bsparam_annotations-1]*MLParamsStruct.bsparam_annolen; % same as above with BS
MLParamsStruct.bsparam_imaxu   = MLParamsStruct.bsparam_imasses + MLParamsStruct.bsparam_masses; % same as above with BS


% a technical constant - SHOULD NOT BE HERE !!!
MLParamsStruct.minimal_log10_t = -10;


end

