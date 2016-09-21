  
if 1
    addpath(genpath('/Users/davidmurphy/GoogleDrive/linked_selection/lsm/cluster_mirror/inf_scripts'))
    label_out_pref = 'dmel_EE2016_';
else
    % cluster mode
    addpath(genpath('/ifs/data/c2b2/gs_lab/dam2214/inf/inf_scripts'))
end

%% load default params and configuration
global MLParamsStruct;
gl_initMLParamsStruct();
fileNames_cluster;
% fileNames_d_mel;
cfg_inf  = LS_DefaultConfiguration(infcfg_file);

%% overwrite config settings where needed (from the fileNames file)

% chromosomes
cfg_inf.chromosomes = 1:length(chr_id);
cfg_inf.inf.chromosomes = cfg_inf.chromosomes;
cfg_inf.chr_id = chr_id;
cfg_inf.inf.chr_id = chr_id;
cfg_inf.output.chromosomes = cfg_inf.chromosomes;
% constants, etc.
cfg_inf.GEs.CalcSW.Ne0 = effective_pop_size;
cfg_inf.GEs.CalcBS.u_del = BS_u_del;
cfg_inf.inf.u_del_max = max_u_del;
cfg_inf.GEs.CalcBS.FE_grid = t_vals;
cfg_inf.GEs.CalcSW.FE_grid = s_vals;
cfg_inf.inf.SWanno2param_mapping = SW_anno_mapping;
cfg_inf.inf.BSanno2param_mapping = BS_anno_mapping;
files_buildGE.outdir = [base_dir GE_dir];
% unfix params for all anno/coeff combinations
for i=1:length(BS_free_params)
  cfg_inf.inf.fixed_params([MLParamsStruct.bsparam_imasses(i) + BS_free_params{i}]) = 0;
end
for i=1:length(SW_free_params)
  cfg_inf.inf.fixed_params([MLParamsStruct.swparam_imasses(i) + SW_free_params{i}]) = 0;
end

%% save structs to files
struct2file(files_buildGE, files_buildGE_file);
struct2file(cfg_inf, infcfg_file);
struct2file(files_invar, files_invar_file);

%% run inference
[calco, mapo, outputo] = LS_InferModel_f( label_out_pref, files_invar_file, files_buildGE_file, files_masks_file, infcfg_file );




