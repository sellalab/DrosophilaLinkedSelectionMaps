
function GEs                     = LS_LoadGridElements( GEs_files, cfg_inf )




skip_generate_BS_maps = 1; 

for c=1:length(cfg_inf.chromosomes)
  for b=1:length(anno_tokens)
    cfgBS_tmplt{c}.cons_table   = Annots{b}{c}.file;
    cfgBS_tmplt{c}.name         = [cfg_inf.genmap_name anno_tokens{b}];
    
    [~, GEs.BSbase{c,b}]  = generateBbase(cfgBS_tmplt{c}, cfg_inf.BSbase_vT, skip_generate_BS_maps);
    
    GEs.BSbase{c,b}.anno           = Annots{b}{c};
    GEs.BSbase{c,b}.chr_id         = Annots{b}{c};
  end
end




end

