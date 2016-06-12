
function calc  = LS_InferModel( outfile_pref, fdata, GEs, cfg_inf, masks, invar_files )

% This function wraps the inference procedure, given prepared input
% structs of variations data (fdata), including grid elements (GEs), an
% inference configuration struct, and if supplied, masks pointing which codons to use. 


global MLParamsStruct;


opts_Alg{1} = optimset('Algorithm','active-set',     'Display','iter');
opts_Alg{2} = optimset('Algorithm','interior-point', 'Display','iter');
opts_Alg{3} = optimset('Algorithm','sqp',            'Display','iter');


%% inits

if isempty( cfg_inf ) % shouldn't be empty...
  te      = LS_DefaultConfiguration( [outfile_pref 'infcfg.txt'], cfg_inf );
  cfg_inf = te.inf;
end



%% set param bounds and set initial values



% genome-wide average basic data statistics, e.g. mean heterozygosity, probability of a segregating site, synonymous divergence.
% NOTE: the statistics are calculated across the whole genome and not just
% for the specific codons used for the inference. I think for our purposes
% here it's good enough, but this may have to be changed.

% fdata = nvdata from LS_InferModel:
gwEstats = genomewideStatisticsLight( fdata(cfg_inf.chromosomes) ); % contains avg poly level, avg het, avg mut rate
% gwEstats = genomewideStatisticsLight_memMap( invar_files, fdata(cfg_inf.chromosomes) ); % MEMORY MAP edition.


% Initializing parameters, either based on data assuming no selection, or from external parameters

tau0 = (gwEstats.Het/gwEstats.MutProx)^-1; % tau0 is the observed average in the empirical data

bnd_t_dist.fixed = cfg_inf.fixed_params; % fixed params set in LS_DefaultConfiguration will carry over here...

% these params are boundaries (H and L mean high and low) for different param values

bnd_t_dist.rangeL(1:MLParamsStruct.length) = 0; % initialize a vector of length(MLParamsStruct) for boundaries
bnd_t_dist.rangeL(MLParamsStruct.tau_pos)  = tau0 / 10; % lower bound for tau is 1/10 observed value
bnd_t_dist.rangeL(MLParamsStruct.swparam_offset+[1:MLParamsStruct.swparam_annotations*(MLParamsStruct.swparam_masses+1)]) = 0; % the lower bound on all SW params set to 0
bnd_t_dist.rangeL(MLParamsStruct.bsparam_offset+[1:MLParamsStruct.bsparam_annotations*(MLParamsStruct.bsparam_masses+1)]) = MLParamsStruct.minimal_log10_t; % ther lower bound on all BS params set to -10
bnd_t_dist.rangeL(MLParamsStruct.bsparam_imaxu) = log10(cfg_inf.u_del_max); % for bsparams, last position in each anno subset min to log u_del_max
bnd_t_dist.rangeL(MLParamsStruct.swparam_imaxu) = 1; % last position in each SW anno subset set to 1 (what is this position for?)

bnd_t_dist.rangeH(MLParamsStruct.tau_pos)  = tau0 * 10; % upper bound on tau set to 10x the emprical tau
bnd_t_dist.rangeH(MLParamsStruct.swparam_offset+[1:MLParamsStruct.swparam_annotations*(MLParamsStruct.swparam_masses+1)]) = 1; % upper bound on all swparams set to 1
bnd_t_dist.rangeH(MLParamsStruct.bsparam_offset+[1:MLParamsStruct.bsparam_annotations*(MLParamsStruct.bsparam_masses+1)]) = log10(10^-6); % I know its -6, just to make things clear % upper bound on bsparams set to -6
bnd_t_dist.rangeH(MLParamsStruct.bsparam_imaxu) = log10(cfg_inf.u_del_max); % set to same as min
bnd_t_dist.rangeH(MLParamsStruct.swparam_imaxu) = 1; % set to same as min...

if isempty(cfg_inf.init_params) % this should be empty under current settings...
  message = 'USING DEFAULT INITIAL PARAMETERS'  
  bnd_t_dist.init                         = bnd_t_dist.rangeL; % set everything to lowest settings as the initial parameters
  bnd_t_dist.init(MLParamsStruct.tau_pos) = tau0;
else
  message = 'READING PRE-DEFINED INITIAL PARAMETERS'
  % bnd_t_dist.init = min( max( cfg_inf.init_params, bnd_t_dist.rangeL ), bnd_t_dist.rangeH );  % edited 01/19/16
  bnd_t_dist.init = min( cfg_inf.init_params, bnd_t_dist.rangeH );

end


% for non-fixed selection weight params, make sure they are not null, since initial null values sometimes stuck the optimization
for k=1:MLParamsStruct.bsparam_annotations
  pidx = intersect( MLParamsStruct.bsparam_imasses(k)+[0:MLParamsStruct.bsparam_masses-1], find(~bnd_t_dist.fixed) ); % find all non-fixed values and set to some starting value greater than 0
  bnd_t_dist.init(pidx) = bnd_t_dist.rangeH;  % EDIT 06/12/16 -- for the non-fixed params, start at the upper bound
end


for k=1:MLParamsStruct.swparam_annotations 
  pidx = intersect( MLParamsStruct.swparam_imasses(k)+[0:MLParamsStruct.swparam_masses-1], find(~bnd_t_dist.fixed) ); % with sw fixed this should return 0 length pidx
  bnd_t_dist.init(pidx) = max( 10^-8, bnd_t_dist.init(pidx) );
end


% fix & null all weights that do not correspond to a vector in GEs (in case we have less than 11 point masses per some annotation)
for k=1:length(cfg_inf.BSanno2param_mapping) % DM edit: change this to run through the length of BSbase instead of bsparam_annotations since BSbase will reflect the final number of annos
  pidx = MLParamsStruct.bsparam_imasses(k)+[length(GEs.BSbase{1,k}.cfg.FE_grid):MLParamsStruct.bsparam_masses-1]; % fix everything after the last grid value in BSbase
  bnd_t_dist.fixed(pidx) = 1;
  bnd_t_dist.init(pidx)  = MLParamsStruct.minimal_log10_t;
end

for k=1:length(cfg_inf.SWanno2param_mapping) % DM edit: this edit is the same as for BSbase above, should be safer this way

  pidx = MLParamsStruct.swparam_imasses(k)+[length(GEs.SWbase{1,k}.cfg.FE_grid):MLParamsStruct.swparam_masses-1]; % fix everything after the last grid value in SWbase
  bnd_t_dist.fixed(pidx) = 1;
  bnd_t_dist.init(pidx)  = 0;
end

%% 
 
% % These are 'just-in-case' old lines, to remind me past mistakes...
%   cfg_inf.rec_th_L = anals_obs{r}.config.rec_th_L;
%   cfg_inf.rec_pol_th_L = cfg_inf.rec_th_L;


  calc         = inferNFSW( GEs, fdata, cfg_inf, bnd_t_dist, masks, [], invar_files);

  struct2file( calc, [outfile_pref 'vparams.txt'] );

  bigo = 7;

end

