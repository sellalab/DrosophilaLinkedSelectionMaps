
function [calc, maps, outcalc]  = LS_InferModel_f( outfile_pref, files_invar_file, files_buildGE_file, files_codonmask_file, cfg_file )

%usage: LS_InferModel_f( 'LS_BS4SW3_', files_invar_file, files_buildGE_file, files_masks_file, infcfg_file );


% inputs:
% outfile_pref          - prefix of inference results output file
% files_invar_file      - configuration file holding file names of input variation data (polymorphism data and mutation rate estimates)
% files_buildGE_file    - configuration file holding file names of pre-calculated grid elements
% files_codonmask_file  - configuration file holding file names of inference mask (i.e. on which codons to run the inference)
% cfg_file              - inference configuration file
% outputs:
% calc                  - a struct holding all of the inference results and related summaries
% maps                  - a struct holding maps of the effects of linked selection 
% outcalc               - a "configuration" struct holding inference results in a way
%                         suitable for saving in a file (saved to outfile_prefXXX file)


% files_invar   = file2struct( files_invar_file );
% files_buildGE = file2struct( files_buildGE_file );
cfg_inf       = file2struct( cfg_file );
cfg_inf.GEs.CalcSW.skip_generate_maps = 1; % temporary(?), just for safety
cfg_inf.GEs.CalcBS.skip_generate_maps = 1; % temporary(?), just for safety


nvdata            = LS_LoadVariationData(    files_invar_file, cfg_inf );
% nvdata            = LS_LoadVariationData_memmap(    files_invar_file, cfg_inf );

GEs               = LS_PrecalcGridElements(  files_buildGE_file, cfg_inf );
% masks             = LS_SetGenomicMask(       files_codonmask_file, cfg_inf );
% calc              = LS_InferModel(           outfile_pref, nvdata, GEs, cfg_inf.inf, masks.inference );
masks             = [];  %here we just provide an empty vector for masks so that the function has enough args:
calc              = LS_InferModel(           outfile_pref, nvdata, GEs, cfg_inf.inf, masks, files_invar_file );  % run without masks, waste of time and space presently...

ii = calc.best_iter;
if ii==-1
  ii = size(calc.params,1);
end
params = calc.params(ii,:);

maps = LS_DrawMap(              outfile_pref, params, GEs, [], cfg_inf );

sparams   = collectParams( calc, GEs, cfg_inf.inf );

outcalc.vparams         = calc.params;
outcalc.params          = sparams.params;
outcalc.fit             = sparams.fit;
% outcalc.stats   = LS_BasicStats( calc, nvdata, GEs, cfg_inf );
outcalc.cfg             = cfg_inf;
outcalc.infiles.var     = file2struct( files_invar_file );
outcalc.infiles.GEs     = file2struct( files_buildGE_file );

% TURNED OFF FOR CLUSTER RUN, CAUSING ERRORS AND NOT NEEDED ANYWAY ... 
% outcalc.infiles.masks   = file2struct( files_codonmask_file );  

struct2file( outcalc, [outfile_pref 'finalresults.txt'] );


end





