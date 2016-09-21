% root dir
if 0
  % local mode
  base_dir = '/Users/davidmurphy/GoogleDrive/linked_selection/lsm/cluster_mirror/';
else
  % cluster mode
  base_dir = '/ifs/data/c2b2/gs_lab/dam2214/inf/';
end

%%
% EVERYTHING BELOW THIS POINT SHOULD WORK IN BOTH CLUSTER MODE AND LOCAL MODE
%%

%% constans for dmel vs. human data
ne_human = 12500.0;  %(based on YRI)
BS_u_del_human = 7.4e-8;
max_udel_human = 1e-7;
ne_dmel    = 1140800.0;
BS_u_del_dmel = 3.5*10^-9;
max_u_del_dmel = 6.8e-9;

%% top level variables
% -> these set additional variables by string combinations, length of arrays, etc

% annotations for bs and cs
BS_anno_tokens = {'primate_cons_exonSegments', 'primate_cons_nonexonSegments'};
% BS_anno_tokens = {'primate_cons95_Segments'};
SW_anno_tokens = {'primate_cons95'};
% there must be at least one annotation for the run to start
assert(~isempty(BS_anno_tokens) | ~isempty(SW_anno_tokens));
% bs and cs coefficients (set to the same range)
t_vals = 10.^-[2:0.5:4.5];
s_vals = 0.01;
% constants (set in defaults, may be reset here)
effective_pop_size = ne_human;
BS_u_del = BS_u_del_human;
max_u_del = max_udel_human;
% annotations directories
BS_anno_dir = 'coords/nr/cons/primate/';
SW_anno_dir = 'subs/cons/';
% pre-calculated grid dirs
GE_dir = 'new_bsmaps';
% genetic map id
genmap_token = 'AA_Map_10kb';
% neutral sites data set to use
% neut_suffix = 'downSample_15pct_neutral';
neut_suffix = 'downSample_1pct_neutral'; % local downsampled version
% mutation rate correction data set
mut_files = 'mutrate';
% chromosome lengths file
chr_features_file = [base_dir 'ch_features/chr_len_all.txt'];

%% dependent param settings

% param mappings per anno (matches number of annotations per bs/cs)
BS_anno_mapping = 1:length(BS_anno_tokens);
SW_anno_mapping = 1:length(SW_anno_tokens);
% free params within the vector of bs and cs params
for i=1:length(BS_anno_mapping)
  BS_free_params{i} = 0:length(t_vals);  % all coeffs for all tokens freed
end
for i=1:length(SW_anno_mapping)
%   SW_free_params{i} = 0:length(s_vals);  % all coeffs for all tokens freed
  SW_free_params{i} = [];  % no params freed
end


%% file paths

% chromosome data
[chr_id, chr_len] = textread(chr_features_file, '%s %d', 'headerlines', 1);
C = length(chr_id);  % num chroms

% genetic map files
for c=1:C
	genmap_files{c} = [base_dir sprintf('maps/%s/%s_%s_window_hg19.txt',genmap_token, chr_id{c}, genmap_token)];
end
% polymorphism data
for c=1:C
  poly_proc_files{c}       = [base_dir sprintf('neut_poly/%s_%s.txt', chr_id{c}, neut_suffix)];
  MutProx_proc_files{c}    = [base_dir sprintf('subs/mutprox/%s_%s.txt', chr_id{c}, mut_files)];
end

% SW_pos_grid_files: these are lists of positions at which to calculate the efects of CS
% ** Note: if neutral sites are filtered for likelihood calculations (e.g. for low cmmb),
% the full set of neutral sites should still be used for calculating the effects of CS
% otherwise SW_pos_grid_files can just be set to the set of neutral sites
for c=1:C
% use positions from neutral poly file:
  SWbase_pos_grid_files{c} = poly_proc_files{c};
  % SWbase_pos_grid_files{c} = [base_dir sprintf('neut_poly/%s_%s.txt', chr_id{c}, neut_suffix)];
end

% annotation files
BS_anno_fileprefs = [base_dir BS_anno_dir];
SW_anno_fileprefs = [base_dir SW_anno_dir];

% default empty cells for each anno
BS_anno1_files = cell(1,C);
BS_anno2_files = cell(1,C);
BS_anno3_files = cell(1,C);
BS_anno4_files = cell(1,C);
BS_anno5_files = cell(1,C);

SW_anno1_files = cell(1,C);
SW_anno2_files = cell(1,C);
SW_anno3_files = cell(1,C);
SW_anno4_files = cell(1,C);

for c=1:C

    % store bs file paths
    if ~isempty(BS_anno_tokens)
      BS_anno1_files{c} = [BS_anno_fileprefs sprintf('ex/%s_%s.bed', chr_id{c}, BS_anno_tokens{1})];
    end
    if length(BS_anno_tokens) > 1
      BS_anno2_files{c} = [BS_anno_fileprefs sprintf('nex/%s_%s.bed', chr_id{c}, BS_anno_tokens{2})];
    end
    if length(BS_anno_tokens) > 2
      BS_anno3_files{c} = [BS_anno_fileprefs sprintf('segs/%s_%s.bed', chr_id{c}, BS_anno_tokens{3})];
    end
    if length(BS_anno_tokens) > 3
      BS_anno4_files{c} = [BS_anno_fileprefs sprintf('segs/%s_%s.bed', chr_id{c}, BS_anno_tokens{4})];
    end
    if length(BS_anno_tokens) > 4
      BS_anno5_files{c} = [BS_anno_fileprefs sprintf('segs/%s_%s.bed', chr_id{c}, BS_anno_tokens{5})];
    end
    % store cs file paths
    if ~isempty(SW_anno_tokens)
      SW_anno1_files{c} = [SW_anno_fileprefs chr_id{c} sprintf('_%s_substitutions.txt', SW_anno_tokens{1})];
    end
    if length(SW_anno_tokens) > 1
      SW_anno2_files{c} = [SW_anno_fileprefs chr_id{c} sprintf('_%s_substitutions.txt', SW_anno_tokens{2})];
    end
    if length(SW_anno_tokens) > 2
      SW_anno3_files{c} = [SW_anno_fileprefs chr_id{c} sprintf('_%s_substitutions.txt', SW_anno_tokens{3})];
    end
    if length(SW_anno_tokens) > 3
      SW_anno4_files{c} = [SW_anno_fileprefs chr_id{c} sprintf('_%s_substitutions.txt', SW_anno_tokens{4})];
    end
end

%% save to config files

infcfg_file = [base_dir sprintf('configs/%scfg.txt', label_out_pref)];
files_invar.poly           = poly_proc_files;
files_invar.mutprox        = MutProx_proc_files;
files_invar.genmap         = genmap_files;
files_invar.genmap_token   = genmap_token;
files_invar.chr_features_file = chr_features_file;
files_invar_file           = [base_dir sprintf('configs/%sfiles_invar.txt', label_out_pref)];
files_buildGE.outfile_pref      = '';
files_buildGE.chr_features_file = chr_features_file;
files_buildGE.chr_id            = chr_id;
files_buildGE.pos_grid_files    = SWbase_pos_grid_files;
files_buildGE.genmap_files      = genmap_files;
files_buildGE.genmap_token      = genmap_token;
% save cs anno file paths
files_buildGE.SW_anno1_files = SW_anno1_files;
files_buildGE.SW_anno2_files = SW_anno2_files;
files_buildGE.SW_anno3_files = SW_anno3_files;
files_buildGE.SW_anno4_files = SW_anno4_files;
files_buildGE.SW_anno_tokens = SW_anno_tokens;
% save bs file paths
files_buildGE.BS_anno1_files = BS_anno1_files;
files_buildGE.BS_anno2_files = BS_anno2_files;
files_buildGE.BS_anno3_files = BS_anno3_files;
files_buildGE.BS_anno4_files = BS_anno4_files;
files_buildGE.BS_anno5_files = BS_anno5_files;
files_buildGE.BS_anno_tokens = BS_anno_tokens;
files_buildGE_file           = [base_dir sprintf('configs/%sfiles_buildGE.txt',label_out_pref)];
files_masks_file             = [base_dir sprintf('configs/%sfiles_masks.txt',label_out_pref)];
% files_bootstrap.chr_features_file = chr_features_file;