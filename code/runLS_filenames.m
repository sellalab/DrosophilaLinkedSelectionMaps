


%% chromosomes

chr_features_file = 'E:/downloads/mel/chromosomal_features_5.25.txt'; %'E:/downloads/mel/chromosomal_features_5.25_autosomes.txt'
[chr_id, chr_len] = textread(chr_features_file, '%s %d', 'headerlines', 1);

C = length(chr_id);

%% genetic maps

genmap_token     = 'Comeron';
genmap_token_alt = 'Song';
for c=1:C
  genmap_files{c}     =  sprintf('E:/downloads/mel/genetic_maps/Comeron/genmap_Comeron_%s.txt',chr_id{c});
  genmap_files_alt{c} =  sprintf('E:/downloads/mel/genetic_maps/Song/RAL_chr%s_Song.txt',chr_id{c});
end



%% selection annotations and LS maps

iSWAnnoExonicNS = 1;
iSWAnnoUTR = 2;
iSWAnnoIntronic = 3;
iSWAnnoIntergenic = 4;

SW_anno_tokens    = {...
  'exonicNS',...
  'UTR',...
  'intronic',...
  'intergenic'};
SW_anno_fileprefs = {...
  'E:/projects/sweeps/NeutFoc/annotations/substitutions_exonicNS_',...
  'E:/projects/sweeps/NeutFoc/annotations/substitutions_UTR_',...
  'E:/projects/sweeps/NeutFoc/annotations/substitutions_intronic_',...
  'E:/projects/sweeps/NeutFoc/annotations/substitutions_intergenic_'};

iBSAnnoExonic = 1;
iBSAnnoUTR = 2;
iBSAnnoIntronic = 3;
iBSAnnoIntergenic = 4;

BS_anno_tokens    = {...
  'exons',...
  'UTRs',...
  'longintrons',...
  'intergenic'};
BS_anno_fileprefs = {...
  'E:/projects/sweeps/NeutFoc/annotations/longest_transcript_533_exon_',...
  'E:/projects/sweeps/NeutFoc/annotations/longest_transcript_533_UTR_',...
  'E:/projects/sweeps/NeutFoc/annotations/longest_transcript_533_longintron_',...
  'E:/projects/sweeps/NeutFoc/annotations/longest_transcript_533_intergenic_'};


for c=1:C
%   for b=1:length(BS_anno_tokens)
    BS_anno1_files{c} = [BS_anno_fileprefs{1} chr_id{c} '.coords'];
    BS_anno2_files{c} = [BS_anno_fileprefs{2} chr_id{c} '.coords'];
    BS_anno3_files{c} = [BS_anno_fileprefs{3} chr_id{c} '.coords'];
    BS_anno4_files{c} = [BS_anno_fileprefs{4} chr_id{c} '.coords'];
%   end
%   for b=1:length(SW_anno_tokens)
    SW_anno1_files{c} = sprintf('E:/projects/sweeps/NeutFoc/annotations/substitutions_%s_%s.txt', chr_id{c}, SW_anno_tokens{1}); %[SW_anno_fileprefs{b} chr_id{c} '.txt'];
    SW_anno2_files{c} = sprintf('E:/projects/sweeps/NeutFoc/annotations/substitutions_%s_%s.txt', chr_id{c}, SW_anno_tokens{2}); %[SW_anno_fileprefs{b} chr_id{c} '.txt'];
    SW_anno3_files{c} = sprintf('E:/projects/sweeps/NeutFoc/annotations/substitutions_%s_%s.txt', chr_id{c}, SW_anno_tokens{3}); %[SW_anno_fileprefs{b} chr_id{c} '.txt'];
    SW_anno4_files{c} = sprintf('E:/projects/sweeps/NeutFoc/annotations/substitutions_%s_%s.txt', chr_id{c}, SW_anno_tokens{4}); %[SW_anno_fileprefs{b} chr_id{c} '.txt'];
    %     SW_anno_files{c}{b} = [SW_anno_fileprefs{b} chr_id{c} '.txt'];
%   end
end


for c=1:C
  SWbase_pos_grid_files{c} = ['E:/projects/sweeps/NeutFoc/annotations/SWbase_positions_' chr_id{c} '.txt'];
end



%% variation data


for c=1:C
  poly_files{c} = {sprintf('E:/downloads/mel/Bailor/strains/Bailor_%s_output1.txt',chr_id{c}),...
    sprintf('E:/downloads/mel/Bailor/strains/Bailor_%s_output2.txt',chr_id{c}),...
    sprintf('E:/downloads/mel/Bailor/strains/Bailor_%s_output3.txt',chr_id{c})};
end

for c=1:C
  poly_proc_files{c}       = sprintf('E:/projects/sweeps/NeutFoc/neut_variation/Filtered_140521_poly_%s.txt',chr_id{c});
  MutProx_proc_files{c}    = sprintf('E:/projects/sweeps/NeutFoc/neut_variation/Filtered_140521_mutrate_%s.txt',chr_id{c});
end

% varStandardPref = 'E:/projects/sweeps/NeutFoc/variation/cleaned';


for c=1:C
  inference_mask_files{c}     = sprintf('E:/projects/sweeps/NeutFoc/pipeline/inference/codons_mask_inference_%s.txt', chr_id{c});
  evaluation_mask_files{c}    = sprintf('E:/projects/sweeps/NeutFoc/pipeline/inference/codons_mask_evaluation_%s.txt', chr_id{c});
  
  bootstrap_output_mask_pref{c}     = sprintf('E:/projects/sweeps/NeutFoc/pipeline/inference/bootstrap/codons_mask_inference_%s', chr_id{c});
end




infcfg_file = 'E:/projects/sweeps/NeutFoc/pipeline/infcfg.txt';

files_invar.poly           = poly_proc_files(1:4);
files_invar.mutprox        = MutProx_proc_files(1:4);
files_invar.genmap         = genmap_files(1:4);
files_invar.genmap_token   = genmap_token;
files_invar_file           = 'E:/projects/sweeps/NeutFoc/pipeline/files_invar.txt';

files_buildGE.outdir            = 'E:/projects/sweeps/NeutFoc/simulated/GEs/'; %'E:/projects/sweeps/NeutFoc/pipeline/GEs/';
files_buildGE.outfile_pref      = '';
files_buildGE.chr_features_file = chr_features_file;
files_buildGE.chr_id            = chr_id(1:4);
files_buildGE.pos_grid_files    = SWbase_pos_grid_files(1:4);
files_buildGE.genmap_files      = genmap_files(1:4);
files_buildGE.genmap_token      = genmap_token;
files_buildGE.SW_anno1_files = SW_anno1_files(1:4);
files_buildGE.SW_anno2_files = SW_anno2_files(1:4);
files_buildGE.SW_anno3_files = SW_anno3_files(1:4);
files_buildGE.SW_anno4_files = SW_anno4_files(1:4);
files_buildGE.SW_anno_tokens = SW_anno_tokens;
files_buildGE.BS_anno1_files = BS_anno1_files(1:4);
files_buildGE.BS_anno2_files = BS_anno2_files(1:4);
files_buildGE.BS_anno3_files = BS_anno3_files(1:4);
files_buildGE.BS_anno4_files = BS_anno4_files(1:4);
files_buildGE.BS_anno_tokens = BS_anno_tokens;
files_buildGE_file           = 'E:/projects/sweeps/NeutFoc/pipeline/files_buildGE.txt';


files_masks.inference           = inference_mask_files(1:4);
files_masks.evaluation          = evaluation_mask_files(1:4);
files_masks_file           = 'E:/projects/sweeps/NeutFoc/pipeline/files_masks.txt';


files_bootstrap.chr_features_file = chr_features_file;
files_bootstrap.input_mask        = files_masks.inference;
files_bootstrap.output_mask_pref  = bootstrap_output_mask_pref(1:4);
files_bootstrap_file                = 'E:/projects/sweeps/NeutFoc/pipeline/inference/bootstrap/files_bootstrap_masks.txt';



