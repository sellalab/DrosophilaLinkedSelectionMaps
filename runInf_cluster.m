
addpath(genpath('/ifs/data/c2b2/gs_lab/dam2214/inf/inf_scripts'))

global MLParamsStruct;
gl_initMLParamsStruct();

%% load basic information: chromosomes data, genetic maps
genmap_id   = 'deCODE';
neut_suffix = '_downSample=25pct_cmmb>0_neutral2.txt';
mut_files   = '_hom2mac_div.txt';

fileNames_cluster

label_out_pref = sprintf('result%s_%s_', neut_suffix, genmap_token);
num_chroms = length(genmap_files);


% this function is unique to each genome, always check the defaults
cfg_inf  = LS_DefaultConfiguration( infcfg_file );

% here the elements of cfg_inf are narrowed down to those from chr1 only:
cfg_inf.inf.SWanno2param_mapping  = [1 2]; % DM edit change the ranges on both of these vectors
cfg_inf.inf.BSanno2param_mapping  = [1 2]; % DM edit change the ranges on both of these vectors 
cfg_inf.inf.chromosomes           = cfg_inf.inf.chromosomes(1:num_chroms); %(chr1)
% cfg_inf.inf.predef_idx_train      = cfg_inf.inf.predef_idx_train(1:num_chroms);
cfg_inf.chr_id                    = cfg_inf.chr_id(1:num_chroms); % change it just for now, easier this way since there is only 1 and it will correctly name files this way
% cfg_inf.predef_idx_test           = cfg_inf.predef_idx_test(1:num_chroms);
cfg_inf.chromosomes               = cfg_inf.chromosomes(1:num_chroms);
cfg_inf.output.chromosomes        = cfg_inf.output.chromosomes(1:num_chroms);

cfg_inf.GEs.CalcSW.FE_grid        = 10.^-1; % DM edit
cfg_inf.GEs.CalcBS.FE_grid        = 10.^-[2:4]; % changed from defaults

%% build Grid Elements (GEs) for BS and SW 

struct2file( files_buildGE, files_buildGE_file );

% this should carry into the precalc script... ? DM edit
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

% struct2file( files_masks, files_masks_file );

[calco, mapo, outputo] = LS_InferModel_f( label_out_pref, files_invar_file, files_buildGE_file, files_masks_file, infcfg_file );


