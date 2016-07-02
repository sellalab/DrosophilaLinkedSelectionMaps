


global MLParamsStruct;

gl_initMLParamsStruct();


%% load basic information: chromosomes data, genetic maps


runLS_filenames;


for c=1:length(genmap_files)
  genmap{c}.file = genmap_files{c};
  genmap{c}.name = genmap_token;
  [genmap{c}.pos, genmap{c}.c, genmap{c}.R]  = textread( genmap{c}.file, '%d\t%f\t%f', 'commentstyle', 'shell' );
end



struct2file( files_buildGE, files_buildGE_file );

cfg_inf  = LS_DefaultConfiguration( infcfg_file );


% TESTING PATCH
cfg_inf.GEs.CalcBS.FE_grid = cfg_inf.GEs.CalcBS.FE_grid(1:3);
cfg_inf.GEs.CalcSW.FE_grid = cfg_inf.GEs.CalcSW.FE_grid(1:3);
cfg_inf.GEs.CalcBS.skip_generate_maps = 0;
cfg_inf.GEs.CalcSW.skip_generate_maps = 0;
struct2file( cfg_inf, infcfg_file );
% TESTING PATCH



%% setup and load Grid Elements (GEs) for BS and SW


% % load grid elements 
cfg_inf.GEs.CalcSW.skip_generate_maps = 1;
cfg_inf.GEs.CalcBS.skip_generate_maps = 1;
struct2file( cfg_inf, infcfg_file );

% [GEs, GE_Annots] = ...
      LS_PrecalcGridElements( files_buildGE_file, infcfg_file );
%         'E:\projects\sweeps\NeutFoc\simulated\GEs',...
%         '',...
%         chr_features_file,...
%         chr_id,...
%         [],...
%         genmap_files,...
%         genmap_token,...
%         SW_anno_files,...
%         SW_anno_tokens,...
%         BS_anno_files,...
%         BS_anno_tokens,...
%         cfg_inf );


% % load grid elements for SW on fitness effects inferred from Sattath et al.'s method
% cfg_inf_STT = cfg_inf;
% cfg_inf_STT.CalcSW.FE_grid = 10.^-[2.1165 5.1269];
% cfg_inf_STT.CalcSW.skip_generate_maps = 1;
% 
% GEs_om = ...
%   LS_PrecalcGridElements( ...
%   'E:\projects\sweeps\NeutFoc\simulated\GEs',...
%   '',...
%   chr_features_file,...
%   chr_id,...
%   pos_grid_files,...
%   genmap_files,...
%   genmap_token,...
%   SW_anno_files,...
%   SW_anno_tokens,...
%   {},...
%   {},...
%   cfg_inf_STT );




% %% load variation data (i.e. load processed data into 'fdata' data-struct)
% 

struct2file( files_invar, files_invar_file );

nvdata = LS_LoadVariationData( files_invar_file );



% create masks of sites used for the inference and for its evaluation
for c=1:length(nvdata)
  chr_mask_notips{c}          = ( (0.05*chr_len(c) < nvdata{c}.Poly.pos)     & (nvdata{c}.Poly.pos     < (1-0.05)*chr_len(c)) );
  chr_mask_c075_rec{c}        = spatAverageRecRat( genmap{c}, nvdata{c}.Poly.pos, cfg_inf.rec_spat_window ) >= 0.75; 
  poly_mask_baseline{c}       = chr_mask_notips{c}         & chr_mask_c075_rec{c};
end
for c=1:length(nvdata)
  cfg_inf.inf.predef_idx_train{c}   = poly_mask_baseline{c};
  cfg_inf.predef_idx_test{c}        = poly_mask_baseline{c};
  f = fopen( files_masks.inference{c}, 'wt' );
  for i=1:length(nvdata{c}.Poly.pos)
  fprintf(f, '%d\t%d\n', nvdata{c}.Poly.pos(i), cfg_inf.inf.predef_idx_train{c}(i) );
  end
  f = fclose(f);
  f = fopen( files_masks.evaluation{c}, 'wt' );
  for i=1:length(nvdata{c}.Poly.pos)
  fprintf(f, '%d\t%d\n', nvdata{c}.Poly.pos(i), cfg_inf.predef_idx_test{c}(i) );
  end
  f = fclose(f);
end


struct2file( files_masks, files_masks_file );






cfg_inf.bootstrap.chromosomes = [1:length(nvdata)];
cfg_inf.bootstrap.repeats = 3;
cfg_inf.bootstrap.block_size = 500000;

struct2file( files_bootstrap, files_bootstrap_file );


cfg_inf = LS_CreateBootstrapMasks( files_bootstrap_file, cfg_inf );


% the next step.. a function collecting bootstrap results; run on cluster; analyzeOBSnCALC 

% [masks, cfg_inf] = LS_SetGenomicMask( files_masks_file, cfg_inf );
% struct2file( cfg_inf, infcfg_file );



%% run inference

[calco, mapo, outputo] = LS_InferModel_f( 'teste', files_invar_file, files_buildGE_file, files_masks_file, infcfg_file );
 
% LS_InferModel_f( files_infer_file, cfg_inf );
% % calc_ex_BS4_SW4         = LS_InferModel( output_pref, fdata, GEs, cfg_inf );



[anals_obs{3}, toymodels]   = LS_StatsDataset( cfg_inf,  nvdata, focals, [], 2 );
anals_BS1234_SW123_a5  = analyzeNFSWcalc( 'BGSnRHH', anals_obs{r}, GEs, fdata, focals, calc_BS1234_SW123_a5{r},  mask_analyses, [], [], [], [] );

  
  






