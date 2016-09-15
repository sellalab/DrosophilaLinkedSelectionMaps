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


%% params
% param mappings
BS_anno_mapping = [1 2 3 4];
SW_anno_mapping = [1 2 3];

% selection coefs
t_vals = 10.^-[1.5:1:5.5];
s_vals = 10.^-[1.5:1:5.5];

% freeing single params at a time for reduced run w/ 3 t-vals
BS_free_params = {[0:4], [0:4], [0:4], [0:4]};
SW_free_params = {[0:4], [0:4], [0:4];


%% file tokens
% annotations for BS and CS in dmel
coords_type = {'exon', 'longintron', 'intergenic', 'UTR'};
substitution_type = {'exonicNS', 'intronic', 'UTR'};

% bs annotations directory
BS_anno_dir = 'coords/dmel/';
CS_anno_dir = 'subs/dmel/'

% LS maps dir
GE_dir = 'new_bsmaps';

% genetic map id
genmap_token = 'genmap_Comeron';

% neutral sites data set to use
neut_prefix = 'Filtered_140521_poly_';

% mutation rate correction data set
mut_prefix = 'Filtered_140521_mutrate_';


%% chromosomes
chr_features_file = [base_dir 'ch_features/dmel_chroms.txt'];
[chr_id, chr_len] = textread(chr_features_file, '%s %d', 'headerlines', 1);

% the number of chromosomes
C = length(chr_id);

%% full file paths 
% write a list of genetic map files including full dir:
for c=1:C
	genmap_files{c} = [base_dir sprintf('maps/comeron/%s_%s.txt', genmap_token, chr_id{c})];
end 

%% BS/CS annos
for i=1:length(substitution_type)
    SW_anno_tokens{i} = substitution_type{i};
end

SW_anno_fileprefs = [base_dir CS_anno_dir];


for i=1:length(coords_type)
    BS_anno_tokens{i} = coords_type{i};  % edited 12/17 -- loop through the coords given on command line externally
end

BS_anno_fileprefs = [base_dir BS_anno_dir];

for c=1:C
    
    BS_anno1_files{c} = '';
    BS_anno2_files{c} = '';
    BS_anno3_files{c} = '';
    BS_anno4_files{c} = '';

    SW_anno1_files{c} = '';
    SW_anno2_files{c} = '';
    SW_anno3_files{c} = '';

    if ~isempty(BS_anno_fileprefs), BS_anno1_files{c} = [BS_anno_fileprefs sprintf('longest_transcript_533_%s_%s.coords' coords_type{1}, chr_id{c})]; end 
    if ~isempty(BS_anno_fileprefs), BS_anno2_files{c} = [BS_anno_fileprefs sprintf('longest_transcript_533_%s_%s.coords' coords_type{2}, chr_id{c})]; end 
    if ~isempty(BS_anno_fileprefs), BS_anno3_files{c} = [BS_anno_fileprefs sprintf('longest_transcript_533_%s_%s.coords' coords_type{3}, chr_id{c})]; end 
    if ~isempty(BS_anno_fileprefs), BS_anno4_files{c} = [BS_anno_fileprefs sprintf('longest_transcript_533_%s_%s.coords' coords_type{4}, chr_id{c})]; end 

    if ~isempty(SW_anno_fileprefs), SW_anno1_files{c} = [SW_anno_fileprefs sprintf('substitutions_%s_%s.txt', chr_id{c} , SW_anno_tokens{1})]; end
    if ~isempty(SW_anno_fileprefs), SW_anno2_files{c} = [SW_anno_fileprefs sprintf('substitutions_%s_%s.txt', chr_id{c} , SW_anno_tokens{2})]; end
    if ~isempty(SW_anno_fileprefs), SW_anno3_files{c} = [SW_anno_fileprefs sprintf('substitutions_%s_%s.txt', chr_id{c} , SW_anno_tokens{3})]; end

end

%% variation data
for c=1:C
  SWbase_pos_grid_files{c} = [base_dir sprintf('neut_poly/dmel/%s_%s.txt', neut_prefix, chr_id{c})];
end

for c=1:C
  poly_proc_files{c}       = [base_dir sprintf('neut_poly/dmel/%s_%s.txt', neut_prefix, chr_id{c})];
  MutProx_proc_files{c}    = [base_dir sprintf('subs/mutprox/dmel/%s_%s.txt', mut_prefix, chr_id{c})]; % DM edit
end

%% write output to a config .txt file
infcfg_file = [base_dir sprintf('configs/%scfg.txt', label_out_pref)];

files_invar.poly           = poly_proc_files;
files_invar.mutprox        = MutProx_proc_files;
files_invar.genmap         = genmap_files;
files_invar.genmap_token   = genmap_token;
files_invar_file           = [base_dir sprintf('configs/%sfiles_invar.txt', label_out_pref)];
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
files_buildGE_file           = [base_dir sprintf('configs/%sfiles_buildGE.txt',label_out_pref)];

files_masks_file             = [base_dir sprintf('configs/%sfiles_masks.txt',label_out_pref)];

files_bootstrap.chr_features_file = chr_features_file;