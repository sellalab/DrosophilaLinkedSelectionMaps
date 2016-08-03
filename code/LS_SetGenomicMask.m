

function masks = LS_SetGenomicMask( codonmask_files, cfg_inf )

% This function loads masks marking which codons to use in the
% inference/evaluation.


Ci = length(cfg_inf.inf.chromosomes);
Ce = length(cfg_inf.chromosomes);


for c=1:Ci
  cfg_inf.inf.inference_file{c} = '';
  masks.inference{c} = [];
end
for c=1:Ce
  cfg_inf.evaluation_file{c} = '';
  masks.evaluation{c}      = [];
end


if isempty(codonmask_files)
  return;
end


if isstr(codonmask_files)
  codonmask_files = file2struct( codonmask_files );
end
  

%% general inits

% [ff_chr_id, ff_chr_len] = textread(chr_features_file, '%s\t%d', 'headerlines', 1);


for c=1:Ci
  if isfield(codonmask_files, 'inference')
%     cfg_inf.inf.inference_file{c} = codonmask_files.inference{c};
    [masks.pos_inference{c}, masks.inference{c}] = textread( codonmask_files.inference{c}, '%d %d' );
  end
end
for c=1:Ce
  if isfield(codonmask_files, 'evaluation')
%     cfg_inf.evaluation_file{c} = codonmask_files.evaluation{c};
    [masks.pos_evaluation{c}, masks.evaluation{c}]      = textread( codonmask_files.evaluation{c}, '%d %d' );
  end
end





end

