  
if 0

    addpath(genpath('/Users/davidmurphy/GoogleDrive/run_lsm/inf_scripts'))
    genmap_id = 'AA_Map';
    neut_suffix = '_mini.txt';
    mut_files = '_mutrate.txt';
    t_vals = [10.^-[2:4]; 10.^-[2:4]; 10.^-[2:4]; 10.^-[2:4]];
    label_out_pref = 'hmm_4an_';
    coords_type = {'phastCons_txn_merged','phastCons_enh_merged', 'phastCons_pro_merged', 'phastCons_reg_merged'};
    substitution_type = {'aa_substitutions_nonseg'};
    map_res = '10kb';
    
end

if 1
    % cluster mode
    addpath(genpath('/ifs/data/c2b2/gs_lab/dam2214/inf/inf_scripts/code/'))
end


global MLParamsStruct;
gl_initMLParamsStruct();

%% load basic information: chromosomes data, genetic maps

fileNames_cluster
%fileNames_cluster_singleAnno

num_chroms = length(genmap_files);


cfg_inf  = LS_DefaultConfiguration( infcfg_file );

cfg_inf.inf.SWanno2param_mapping  = 1;  % ONE MAPPING FOR NS SUBSTITUTIONS
cfg_inf.inf.BSanno2param_mapping  = BS_anno_mapping;  % 01/26/16 -- switch off for one bs param w/ phastcons

%% EDIT THE FILES DIR IF NEEDED:
% files_buildGE.outdir = [base_dir 'bsmaps_gcons'];
files_buildGE.outdir            = [base_dir GE_dir];

% THESE LINES ARE NO LONGER NEEDED FOR RUNNING ALL CHROMS... (comment out 01/19/16)
% cfg_inf.inf.chromosomes           = cfg_inf.inf.chromosomes(1:num_chroms); %(chr1)
% cfg_inf.chr_id                    = cfg_inf.chr_id(1:num_chroms); % change it just for now, easier this way since there is only 1 and it will correctly name files this way
% cfg_inf.chromosomes               = cfg_inf.chromosomes(1:num_chroms);
% cfg_inf.output.chromosomes        = cfg_inf.output.chromosomes(1:num_chroms);

% cfg_inf.GEs.CalcSW.FE_grid        = 10.^-3;
% edit 12/17 -- now FE grid has cells for each anno -- must revise formatting to match the following:
%	cfg_inf.GEs.CalcBS.FE_grid{1}     = 10.^-[2:4];
%	cfg_inf.GEs.CalcBS.FE_grid{2}     = 10.^-[2:4];

%% PROVIDE INITIAL PARAMS:

% get them from a file:
% if 0  % switch
% 	params_best = '8pt_retry_best_params.txt';
% 	ff = fopen(params_best, 'r');
% 	best_params = fscanf(ff, '%f', [1, 106]);
% 	cfg_inf.inf.init_params = best_params;
% end

%% BS GRID:
cfg_inf.GEs.CalcBS.FE_grid = t_vals; % define based on t_vals cells


%% SET FIXED PARAMS
for i=1:length(BS_free_params)
  cfg_inf.inf.fixed_params([MLParamsStruct.bsparam_imasses(i) + BS_free_params{i}]) = 0;
  % cfg_inf.inf.fixed_params([MLParamsStruct.bsparam_imasses(1) + BS_free_params1]) = 0;
  % cfg_inf.inf.fixed_params([MLParamsStruct.bsparam_imasses(2) + BS_free_params2]) = 0;
  % cfg_inf.inf.fixed_params([MLParamsStruct.bsparam_imasses(4) + BS_free_params4]) = 0;
end
  % cfg_inf.inf.fixed_params([MLParamsStruct.bsparam_imasses(3) + [0:5]]) = 0;  % 03/03/16 two more sets of params for chromHMM annos
  % cfg_inf.inf.fixed_params([MLParamsStruct.bsparam_imasses(4) + [0:5]]) = 0;



%% build Grid Elements (GEs) for BS and SW 
struct2file( files_buildGE, files_buildGE_file );

cfg_inf.GEs.CalcSW.skip_generate_maps = 1; % DM edit turn on or off here for as the final setting
cfg_inf.GEs.CalcBS.skip_generate_maps = 1; % depends on whether or not they are ready (slow step)

% after this setting the final struct file is created:
struct2file( cfg_inf, infcfg_file );

%files_buildGE_file, infcfg_file, chr_features_file %debug printouts for LS_PrecalcGridElements()

% this can be skipped because 'LS_InferModel_f' builds the GEs if configured to
% LS_PrecalcGridElements( files_buildGE_file, infcfg_file ); % DM edit-- turn off if maps complete



%% run inference

% since we built the GEs, don't do it again - just load them
cfg_inf.GEs.CalcSW.skip_generate_maps = 1;
cfg_inf.GEs.CalcBS.skip_generate_maps = 1;
struct2file( cfg_inf, infcfg_file );

struct2file( files_invar, files_invar_file );


[calco, mapo, outputo] = LS_InferModel_f( label_out_pref, files_invar_file, files_buildGE_file, files_masks_file, infcfg_file );




