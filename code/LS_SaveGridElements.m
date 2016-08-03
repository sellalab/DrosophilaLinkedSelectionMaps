
function LS_SaveGridElements( ourput_pref, GEs, cfg_inf )

% chr_id = {'2L', '2R', '3L', '3R', 'X'};


for c=1:size(GEs.BSbase,1) % chromosomes
  for b=1:size(GEs.BSbase,2) % annotations
    for t=1:length(GEs.BSbase{c,b}.Lj) % annotations
    % write header with relevant parameters
    
    struct2file();
    
    % write BS values table
    f = fopen( , 'at' );
    
    f = fclose(f);
    end
  end
end





skip_generate_BS_maps = 1; 

for c=1:length(cfg_inf.chromosomes)
  for b=1:length(anno_tokens)
    cfgBS_tmplt{c}.cons_table   = Annots{b}{c}.file;
    cfgBS_tmplt{c}.name         = [cfg_inf.genmap_name anno_tokens{b}];
    
    [~, GEs.BSbase{c,b}]  = generateBbase(cfgBS_tmplt{c}, cfg_inf.BSbase_vT, skip_generate_BS_maps);
    
    GEs.BSbase{c,b}.anno           = Annots{b}{c};
  end
end




end

