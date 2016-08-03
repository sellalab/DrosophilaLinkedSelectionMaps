% root dir
if 0
  % local mode
  base_dir = '/Users/davidmurphy/GoogleDrive/linked_selection/lsm/cluster_mirror/';
end

if 1
  % cluster mode
  base_dir = '/ifs/data/c2b2/gs_lab/dam2214/inf/';
end

%%
% EVERYTHING BELOW THIS POINT SHOULD WORK IN BOTH CLUSTER MODE AND LOCAL MODE
%%

% params to change depending on annos and number of annos
% coords_type = {'phastCons_codingSegments', 'phastCons_utrSegments', 'phastCons_nearExonSegments', 'phastCons_nonPeriExonicSegments'}
% coords_type = {'random_1', 'random_2'}
% coords_type = {'primateCons_exonSegments', 'primateCons_nonexonicSegments'}
% coords_type = {'phastCons_exonSegments', 'phastCons_nonexonicSegments'}

% --- exon/nonexon partition
% coords_type = {'ape_cons_exonSegments', 'ape_cons_nonexonSegments'}
% coords_type = {'primate_cons_exonSegments', 'primate_cons_nonexonSegments'}
% coords_type = {'prosimian_cons_exonSegments', 'prosimian_cons_nonexonSegments'}
% coords_type = {'euarchontoglires_cons_exonSegments', 'euarchontoglires_cons_nonexonSegments'}
% coords_type = {'laurasiatheria_cons_exonSegments', 'laurasiatheria_cons_nonexonSegments'}
% coords_type = {'afrotheria_cons_exonSegments', 'afrotheria_cons_nonexonSegments'}
% coords_type = {'mammal_cons_exonSegments', 'mammal_cons_nonexonSegments'}
% coords_type = {'birds_cons_exonSegments', 'birds_cons_nonexonSegments'}
% coords_type = {'fish_cons_exonSegments', 'fish_cons_nonexonSegments'}
% coords_type = {'lamprey_cons_exonSegments', 'lamprey_cons_nonexonSegments'}

% --- no partition
% coords_type = {'ape_cons98_Segments'}
% coords_type = {'primate_cons_exnex_merge'}
coords_type = {'primate_cons95_Segments'}
% coords_type = {'prosimian_cons98_Segments'}
% coords_type = {'euarchontoglires_cons90_Segments'}
% coords_type = {'laurasiatheria_cons90_Segments'}
% coords_type = {'afrotheria_cons90_Segments'}
% coords_type = {'mammal_cons90_Segments'}
% coords_type = {'birds_cons90_Segments'}
% coords_type = {'fish_cons90_Segments'}
% coords_type = {'lamprey_cons90_Segments'}
% coords_type = {'cadd_cons'}


% --- coords from tss, utr5, cds, peri, utr3:
% coords_type = {'knownGene_nonredundant_tssSegments', 'knownGene_nonredundant_utr5Segments', 'knownGene_nonredundant_codingSegments', 'knownGene_nonredundant_periSegments', 'knownGene_nonredundant_utr3Segments'}

% --- chromHMM coords
% coords_type = {'H1hesc_enh_merged', 'H1hesc_pro_merged', 'H1hesc_txn_merged', 'H1hesc_ins_merged', 'H1hesc_rep_merged'}
% coords_type = {'H1hesc_enh_merged', 'H1hesc_pro_merged', 'H1hesc_txn_merged', 'H1hesc_ins_merged'}
% coords_type = {'All_enh_merged', 'All_pro_merged', 'All_txn_merged', 'All_ins_merged'}

substitution_type = {'aa_substitutions_nonseg'}

% --- bs annotations directory
% BS_anno_dir = 'coords/nr/'
% BS_anno_dir = 'coords/nr/cons/ape/'
BS_anno_dir = 'coords/nr/cons/primate/'
% BS_anno_dir = 'coords/nr/cons/prosimian/'
% BS_anno_dir = 'coords/nr/cons/euarchontoglires/'
% BS_anno_dir = 'coords/nr/cons/laurasiatheria/'
% BS_anno_dir = 'coords/nr/cons/afrotheria/'
% BS_anno_dir = 'coords/nr/cons/mammal/'
% BS_anno_dir = 'coords/nr/cons/birds/'
% BS_anno_dir = 'coords/nr/cons/fish/'
% BS_anno_dir = 'coords/nr/cons/lamprey/'
% BS_anno_dir = 'coords/hmm/'
% BS_anno_dir = 'coords/cadd/'


% --- bs maps dir
GE_dir = 'new_bsmaps'

% genetic map id
genmap_id = 'AAmod'
map_res = '2kb'

% neutral sites data set to use
neut_suffix = '_downSample_15pct_neutral.txt'

% mutation rate correction data set
mut_files = '_mutrate.txt'

% BS_anno_mapping = [1 2 3 4 5];
% BS_anno_mapping = [1 2 3 4];
% BS_anno_mapping = [1 2];
BS_anno_mapping = 1;

% --- selection coefs
% t_vals = [10.^-[2:0.5:4.5]; 10.^-[2:0.5:4.5]; 10.^-[2:0.5:4.5]; 10.^-[2:0.5:4.5]; 10.^-[2:0.5:4.5]];
% t_vals = [10.^-[2:0.5:4.5]; 10.^-[2:0.5:4.5]; 10.^-[2:0.5:4.5]; 10.^-[2:0.5:4.5]];
% t_vals = [10.^-[1:0.5:4.5]; 10.^-[1:0.5:4.5]; 10.^-[1:0.5:4.5]; 10.^-[1:0.5:4.5]];
% t_vals = [10.^-[1:4]; 10.^-[1:4]; 10.^-[1:4]; 10.^-[1:4]];
% t_vals = [10.^-[1:0.5:4.5]; 10.^-[1:0.5:4.5]];
% t_vals = [10.^-[2:0.5:4.5]; 10.^-[2:0.5:4.5]];
t_vals = [10.^-[2:0.5:4.5]]
% t_vals = [10.^-[2:0.25:4.5]]
% t_vals = [10.^-[2:4]]

% freeing single params at a time for reduced run w/ 3 t-vals
% BS_free_params = {[0 2]};


% BS_free_params = {[0:2]};
BS_free_params = {[0:5]};
% BS_free_params = {[0:5], [0:5]};
% BS_free_params = {[0:5], [0:5], [0:5], [0:5]};
% BS_free_params = {[0:5], [0:5], [0:5], [0:5], [0:5]};

%% chromosomes
chr_features_file = [base_dir 'ch_features/chr_len_all.txt'];
[chr_id, chr_len] = textread(chr_features_file, '%s %d', 'headerlines', 1);

% the number of chromosomes
C = length(chr_id)

%% genetic maps
genmap_token = sprintf('%s_%s', genmap_id, map_res);

% write a list of genetic map files including full dir:
for c=1:C
	genmap_files{c} = [base_dir sprintf('maps/%s/%s_%s_window_hg19.txt',genmap_id, chr_id{c}, genmap_token)];
end 


%% selection annotation labels:

% modified to be completely shut off 9/28/15
% SW_anno_tokens = {...
% 'utr5_nonsynonymous',... % these correspond to the .sw files, must match
% 'non_coding' };


for i=1:length(substitution_type)
    SW_anno_tokens{i} = substitution_type{i};
end

SW_anno_fileprefs = [base_dir 'subs/derived_substitutions/'];

for i=1:length(coords_type)
    BS_anno_tokens{i} = coords_type{i};  % edited 12/17 -- loop through the coords given on command line externally
end

BS_anno_fileprefs = [base_dir BS_anno_dir];

for c=1:C
    
    BS_anno1_files{c} = '';
    BS_anno2_files{c} = '';
    BS_anno3_files{c} = '';
    BS_anno4_files{c} = '';
    BS_anno5_files{c} = '';


%    SW_anno1_files{c} = '';
%    SW_anno2_files{c} = '';

    % if ~isempty(BS_anno_fileprefs), BS_anno1_files{c} = [BS_anno_fileprefs sprintf('mergecons/%s_%s.bed', chr_id{c}, coords_type{1})]; end 
    if ~isempty(BS_anno_fileprefs), BS_anno1_files{c} = [BS_anno_fileprefs sprintf('segs/%s_%s.bed', chr_id{c}, coords_type{1})]; end 
    % if ~isempty(BS_anno_fileprefs), BS_anno1_files{c} = [BS_anno_fileprefs sprintf('%s_%s.bed', chr_id{c}, coords_type{1})]; end 

    % if ~isempty(BS_anno_fileprefs), BS_anno1_files{c} = [BS_anno_fileprefs sprintf('ex/%s_%s.bed', chr_id{c}, coords_type{1})]; end 
    % if ~isempty(BS_anno_fileprefs), BS_anno2_files{c} = [BS_anno_fileprefs sprintf('nex/%s_%s.bed', chr_id{c}, coords_type{2})]; end

    % if ~isempty(BS_anno_fileprefs), BS_anno1_files{c} = [BS_anno_fileprefs sprintf('enh/%s_%s.bed', chr_id{c}, coords_type{1})]; end 
    % if ~isempty(BS_anno_fileprefs), BS_anno2_files{c} = [BS_anno_fileprefs sprintf('pro/%s_%s.bed', chr_id{c}, coords_type{2})]; end
    % if ~isempty(BS_anno_fileprefs), BS_anno3_files{c} = [BS_anno_fileprefs sprintf('txn/%s_%s.bed', chr_id{c}, coords_type{3})]; end
    % if ~isempty(BS_anno_fileprefs), BS_anno4_files{c} = [BS_anno_fileprefs sprintf('ins/%s_%s.bed', chr_id{c}, coords_type{4})]; end
    % if ~isempty(BS_anno_fileprefs), BS_anno5_files{c} = [BS_anno_fileprefs sprintf('rep/%s_%s.bed', chr_id{c}, coords_type{5})]; end

    % if ~isempty(BS_anno_fileprefs), BS_anno1_files{c} = [BS_anno_fileprefs sprintf('tss/%s_%s.bed', chr_id{c}, coords_type{1})]; end 
    % if ~isempty(BS_anno_fileprefs), BS_anno2_files{c} = [BS_anno_fileprefs sprintf('utr5/%s_%s.bed', chr_id{c}, coords_type{2})]; end
    % if ~isempty(BS_anno_fileprefs), BS_anno3_files{c} = [BS_anno_fileprefs sprintf('cds/%s_%s.bed', chr_id{c}, coords_type{3})]; end
    % if ~isempty(BS_anno_fileprefs), BS_anno4_files{c} = [BS_anno_fileprefs sprintf('peri/%s_%s.bed', chr_id{c}, coords_type{4})]; end
    % if ~isempty(BS_anno_fileprefs), BS_anno5_files{c} = [BS_anno_fileprefs sprintf('utr3/%s_%s.bed', chr_id{c}, coords_type{5})]; end

    if ~isempty(SW_anno_fileprefs), SW_anno1_files{c} = [SW_anno_fileprefs chr_id{c} sprintf('_%s.txt', SW_anno_tokens{1})]; end
%     if ~isempty(SW_anno_fileprefs), SW_anno2_files{c} = [SW_anno_fileprefs chr_id{c} '.subst.nc'];, end

end

% prepare CS grid file names (these are lists of positions at which to calculate the efects of CS, instead of using a fixed-distance grid)
for c=1:C
% use positions from neutral poly file:
  SWbase_pos_grid_files{c} = [base_dir sprintf('neut_poly/%s%s', chr_id{c}, neut_suffix)];
end



%% variation data

for c=1:C
  poly_proc_files{c}       = [base_dir sprintf('neut_poly/%s%s', chr_id{c}, neut_suffix)];
  MutProx_proc_files{c}    = [base_dir sprintf('subs/mutprox/%s%s', chr_id{c}, mut_files)]; % DM edit
end

%% write output to a config .txt file
infcfg_file = [base_dir sprintf('configs/%s_cfg.txt', label_out_pref)];

files_invar.poly           = poly_proc_files;
files_invar.mutprox        = MutProx_proc_files;
files_invar.genmap         = genmap_files;
files_invar.genmap_token   = genmap_token;
files_invar_file           = [base_dir sprintf('configs/%s_files_invar.txt', label_out_pref)];
% files_buildGE.outdir            = [base_dir 'human_GEs']; 
% files_buildGE.outdir            = [base_dir 'bsmaps_gcons'];  % edited 01/15/16 to try old maps again (edit 01/29/16 -- set in main)
files_buildGE.outfile_pref      = '';
files_buildGE.chr_features_file = chr_features_file;
files_buildGE.chr_id            = chr_id; 
files_buildGE.pos_grid_files    = SWbase_pos_grid_files;
files_buildGE.genmap_files      = genmap_files; 
files_buildGE.genmap_token      = genmap_token;
files_buildGE.SW_anno1_files = SW_anno1_files;  % edited 12/17 (turn on)
% files_buildGE.SW_anno2_files = SW_anno2_files; 
files_buildGE.SW_anno_tokens = SW_anno_tokens;  % edited 12/27
files_buildGE.BS_anno1_files = BS_anno1_files;
% files_buildGE.BS_anno2_files = BS_anno2_files;  % edited 12/17
% files_buildGE.BS_anno3_files = BS_anno3_files;  % 03/03/16 adding MORE annotations with chromHMM
% files_buildGE.BS_anno4_files = BS_anno4_files;
% files_buildGE.BS_anno5_files = BS_anno5_files;
files_buildGE.BS_anno_tokens = BS_anno_tokens;
files_buildGE_file           = [base_dir sprintf('configs/%s_files_buildGE.txt',label_out_pref)];

files_masks_file             = [base_dir sprintf('configs/%s_files_masks.txt',label_out_pref)];

files_bootstrap.chr_features_file = chr_features_file;