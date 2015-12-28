
function ncfg   = LS_DefaultConfiguration( outfile, cfg )

% This function creates a new full inference configuration struct or
% updates/completes a given one.

% inputs:
% outfile    - output configuration file name
% cfg        - existing configuration struct
% outputs:
% ncfg       - new/updated configuration struct

% called explicitly in the main script also
global MLParamsStruct;
gl_initMLParamsStruct();

% Initialize all configuration fields (for both data processing, model inferece and model evaluation).
% If 'cfg' has some of the fields, only the missing ones are initialized.

if nargin < 2
  ncfg = [];
else
  ncfg = cfg;
end

if isstr(ncfg)
  ncfg = file2struct(ncfg);
end

if ~isfield(ncfg,     'inf'),                   ncfg.inf= []; ,end
if ~isfield(ncfg, 'predef_idx_test'),            ncfg.predef_idx_test = {[],[],[],[]}; ,end % on which sites to run inference, subsets...
if ~isfield(ncfg, 'chromosomes'),        ncfg.chromosomes = [1:22]; ,end 
if ~isfield(ncfg,        'GEs'),          ncfg.GEs = []; ,end


%% Ask Eyal about these settings:
if ~isfield(ncfg.inf, 'tau_sweeps'),            ncfg.inf.tau_sweeps = -1; ,end % when this is -1, calculate by scaling tau_mutprox, using the parameters 'ratio_tau_mutprox_tau_sweeps'
if ~isfield(ncfg.inf, 'ratio_tau_mutprox_tau_sweeps'), ncfg.inf.ratio_tau_mutprox_tau_sweeps = 3.25; ,end % default = 1 ?  % UNCLEAR WHAT THIS DOES...
if ~isfield(ncfg.inf, 'remote2closeDiv'),       ncfg.inf.remote2closeDiv = 3.2583; ,end % the only current usage of this param is to calculate "tau_subs"

%% inference params
% this is still below mcvicker's inferred u_del for exons:
if ~isfield(ncfg.inf, 'u_del_max'),             ncfg.inf.u_del_max = 6*10^-8; ,end % DM edit: make the maximal deleterious mutation rate MUCH higher...
if ~isfield(ncfg.inf, 'predef_idx_train'),            ncfg.inf.predef_idx_train = {[],[],[],[]}; ,end
if ~isfield(ncfg.inf, 'init_params'),        ncfg.inf.init_params = []; ,end
if ~isfield(ncfg.inf, 'optim_algs'),            ncfg.inf.optim_algs = [1 2 3]; ,end
if ~isfield(ncfg.inf, 'opt_type'), ncfg.inf.opt_type = {'active-set', 'interior-point', 'sqp'};   ,end
if ~isfield(ncfg.inf, 'fmincon_retries '),      ncfg.inf.fmincon_retries = 0; ,end
if ~isfield(ncfg.inf, 'chromosomes'),           ncfg.inf.chromosomes = [1:22]; ,end % DM edit for human chroms
if ~isfield(ncfg.inf, 'constraint_u'),          ncfg.inf.constraint_u = 0; ,end % DM edit, set to 0 (no constraint on u_del
if ~isfield(ncfg.inf, 'use_fake_subs'),         ncfg.inf.use_fake_subs     = 0;,end % DM edit set these to 0
if ~isfield(ncfg.inf, 'SWanno2param_mapping'),  ncfg.inf.SWanno2param_mapping = [1 2]; ,end % DM edit for just coding/noncoding
if ~isfield(ncfg.inf, 'BSanno2param_mapping'),  ncfg.inf.BSanno2param_mapping = [1 2]; ,end % DM edit for just coding/noncoding

%% fixing and freeing different params:
if ~isfield(ncfg.inf, 'fixed_params') 

    % Set to 0 anything that will be allowed to vary in the composite likelihood maximization:
    ncfg.inf.fixed_params(1:MLParamsStruct.length) = 1; % fix all parameters in the initial config
    ncfg.inf.fixed_params(MLParamsStruct.tau_pos ) = 0; % TAU: ratio of 1/(pi0): never fix this or the inference won't work 

    ncfg.inf.fixed_params([MLParamsStruct.bsparam_imasses(1) + [0:2]]) = 0;
    ncfg.inf.fixed_params([MLParamsStruct.bsparam_imasses(2) + [0:2]]) = 0;
    
end

%% string names for chroms:
if ~isfield(ncfg.inf, 'chr_id'),   
    for i=1:22 % DM edit chromosome number
        ncfg.chr_id{i} = sprintf('chr%d', i);
    end
end

%% recombination thresholds --> NOT USED IN CURRENT SETTINGS
if ~isfield(ncfg, 'foc_rec_th_L'),        ncfg.foc_rec_th_L = 0.75; ,end % DM Note: not used anywhere in current implementation. (cM/Mb min choice of focal subs for sweeps )
if ~isfield(ncfg, 'sort_strand_dir'),        ncfg.sort_strand_dir = 0; ,end % DM edit turn off
if ~isfield(ncfg, 'foc_rec_th_H'),        ncfg.foc_rec_th_H = Inf; ,end % upper bound for rec
if ~isfield(ncfg, 'rec_th_L'),        ncfg.rec_th_L = 0.75; ,end
if ~isfield(ncfg, 'rec_th_H'),        ncfg.rec_th_H = Inf; ,end
if ~isfield(ncfg, 'rec_pol_th_L'),        ncfg.rec_pol_th_L = 0.75; ,end
if ~isfield(ncfg, 'rec_pol_th_H'),        ncfg.rec_pol_th_H = Inf; ,end

%% plotting settings
if ~isfield(ncfg, 'gL_COL'),        ncfg.gL_COL = 0.0011; ,end
if ~isfield(ncfg, 'gdelta_COL'),        ncfg.gdelta_COL = 10^-6; ,end
if ~isfield(ncfg, 'halfway'),        ncfg.halfway = 0; ,end
if ~isfield(ncfg, 'collated_weights'),        ncfg.collated_weights = 0; ,end
if ~isfield(ncfg, 'scale_dS'),        ncfg.scale_dS = 0.4; ,end % DM Note: commented out in all inference scripts. dS is inferred it appears
if ~isfield(ncfg, 'use_paml_dS'),        ncfg.use_paml_dS = 1; ,end % DM Note: commented out in downstream scripts
if ~isfield(ncfg, 'min_paml_codons'),        ncfg.min_paml_codons = 50; ,end % DM Note: not used in inference
if ~isfield(ncfg, 'smoothgenmap'),        ncfg.smoothgenmap = 1; ,end
if ~isfield(ncfg, 'smoothgenmap_win'),        ncfg.smoothgenmap_win = 10^-6; ,end
if ~isfield(ncfg, 'rec_spat_window'),        ncfg.rec_spat_window = 10^5; ,end % DM Note: turned off in current implementation

%% settings for the output
% if ~isfield(ncfg, 'genmap_name'),        ncfg.genmap_name = 'AA_Map'; ,end % REMOVE THIS SECTION, NAME IS DEFINED ELSEWHERE
if ~isfield(ncfg, 'output'),             ncfg.output = []; ,end
if ~isfield(ncfg.output, 'LSmap_res'),          ncfg.output.LSmap_res = 100; ,end
if ~isfield(ncfg.output, 'chromosomes'),          ncfg.output.chromosomes = [1:22]; ,end % DM edit chrom number
if ~isfield(ncfg.output, 'spatial_resolution'),          ncfg.output.spatial_resolution = 1000; ,end


%% Sweep GE inits:
if ~isfield(ncfg.GEs,        'CalcSW'),          ncfg.GEs.CalcSW = []; ,end
if ~isfield(ncfg.GEs.CalcSW, 'StopSum'),         ncfg.GEs.CalcSW.StopSum       = 1;              ,end
if ~isfield(ncfg.GEs.CalcSW, 'InterpMethod'),    ncfg.GEs.CalcSW.InterpMethod  = 'linear';       ,end
if ~isfield(ncfg.GEs.CalcSW, 'trap_aprx'),       ncfg.GEs.CalcSW.trap_aprx     = 'diffusion';    ,end  %DURRETT1;
if ~isfield(ncfg.GEs.CalcSW, 'gMaxDistScaled'),  ncfg.GEs.CalcSW.gMaxDistScaled= 1;              ,end  %scale the radius of maximal summation in proportion to S
if ~isfield(ncfg.GEs.CalcSW, 'gMaxDist'),        ncfg.GEs.CalcSW.gMaxDist      = 1;              ,end  %if no scaling, absolute distance in morgans, otherwise the proportionality constant between S and maximal radius (0.1 is reasonable given the rule-of-thumb r=0.1s, 1 is very safe)
if ~isfield(ncfg.GEs.CalcSW, 'FE_grid'),         ncfg.GEs.CalcSW.FE_grid = 10.^-[1:0.25:3.5];    ,end % DM edit minimum set to 10^-3.5 for more realistic assumptions in humans
if ~isfield(ncfg.GEs.CalcSW, 'Ne0'),             ncfg.GEs.CalcSW.Ne0 = 1.25*10^4;              ,end % DM edit for human Ne0 ~1.25e4
if ~isfield(ncfg.GEs.CalcSW, 'skip_generate_maps'),  ncfg.GEs.CalcSW.skip_generate_maps = 0; ,end % DM edit turn off
if ~isfield(ncfg.GEs.CalcSW, 'min_mapdist_SW_spat_grid'), ncfg.GEs.CalcSW.min_mapdist_SW_spat_grid = 3*10^-9; ,end  % this setting influences the spacing of grid elements, this is likely to be important for the change to human data
 
%% BS GE inits:
if ~isfield(ncfg.GEs,        'CalcBS'),          ncfg.GEs.CalcBS = []; ,end
% if ~isfield(ncfg.GEs.CalcBS, 'FE_grid'),             ncfg.GEs.CalcBS.FE_grid = 10.^-[1:0.25:3.5];          ,end % DM edit min set to 10^-3.5 for more realistic assumptions in humans

% conditions where there are TWO different ranges of FE_grids:
% if ~isfield(ncfg.GEs.CalcBS, 'FE_grid'), ncfg.GEs.CalcBS.FE_grid{1} = 10.^-[2:4]; end  % exonic
% ncfg.GEs.CalcBS.FE_grid{2} = 10.^-[4.5:0.5:5.5];  % nonexonic
% TENTATIVE:
if ~isfield(ncfg.GEs.CalcBS, 'FE_grid'), ncfg.GEs.CalcBS.FE_grid = [10.^-[2.5:0.5:4]; 10.^-[4.5:0.5:6]]; end

if ~isfield(ncfg.GEs.CalcBS, 'u_del'),               ncfg.GEs.CalcBS.u_del      = 7.4*10^-8;           ,end % DM edit: use human mutation rate from Leffler review. This is used in the bkgd program so it has an effect on the results
if ~isfield(ncfg.GEs.CalcBS, 'skip_generate_maps'),  ncfg.GEs.CalcBS.skip_generate_maps = 0; ,end % DM edit turn off
if ~isfield(ncfg.GEs.CalcBS, 'B_res'),               ncfg.GEs.CalcBS.B_res = 100; ,end
if ~isfield(ncfg.GEs.CalcBS, 'exe_file'),            ncfg.GEs.CalcBS.exe_file = '/Users/davidmurphy/GoogleDrive/linked_selection/lsm/cpp/bkgd_map/src/calc_bkgd'; ,end % DM edit
if ~isfield(ncfg.GEs.CalcBS, 't_dist_type'),       ncfg.GEs.CalcBS.t_dist_type = 'POINT'; ,end

if ~isempty(outfile)
  struct2file( ncfg, outfile );
end


end

