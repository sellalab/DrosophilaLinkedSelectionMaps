
function cfg = LS_CreateBootstrapMasks( files_bootstrap, cfg )

if isstr(files_bootstrap)
  files_bootstrap = file2struct( files_bootstrap );
end

[ff_chr_id, ff_chr_len] = textread( files_bootstrap.chr_features_file, '%s %d', 'commentstyle', 'shell' ); %, 'headerlines', 1
for c=1:length(cfg.bootstrap.chromosomes)
  ichr(c) = find(strcmp( ff_chr_id, cfg.chr_id{cfg.bootstrap.chromosomes(c)} ));
  chr_len(c) = ff_chr_len(ichr(c));
  chr_id(c)  = ff_chr_id( ichr(c));
end


% read 0-1 mask of which codons included in the analysis

% if isfield(files_bootstrap, 'inference')
Ci = length(files_bootstrap.input_mask); %files_bootstrap.inference
for c=1:Ci
  [pos_input{c}, m_input{c}] = textread( files_bootstrap.input_mask{c}, '%d\t%d' );
end
% end

% if isfield(files_bootstrap,'evaluation')
%   Ce = length(files_bootstrap.evaluation);
%   for c=1:Ce
%     m_evaluation{c} = textread( files_bootstrap.evaluation{c}, '%d' );
%   end
% end


% divide mask to M equally spaced windows

nblocks   = 0;
block_chr = [];
block_idx = [];
for c=1:Ci
  nblocks_chr = floor( chr_len(c) / cfg.bootstrap.block_size ); % I leave out the last less-than-a-full block
  block_chr  = [block_chr c*ones([1 nblocks_chr])];
  block_idx  = [block_idx [1:nblocks_chr]];
  nblocks = nblocks + nblocks_chr;
end


% initialize RNG
if isfield( cfg.bootstrap, 'seed' ) & cfg.bootstrap.seed>=0
  rng(cfg.bootstrap.seed);
else
  rng('shuffle');
end
rng_state = rng;
cfg.bootstrap.seed = rng_state.Seed;


% mode 1: bootstrap 'repeats' masks, each containing 'nblock' windows with replacement

for i=1:cfg.bootstrap.repeats
  
  % sample nblocks with replacement from the range 1:nblocks (=bootstrapping)
  iblocks = randi(nblocks, 1, nblocks);
  
  % initialize output mask
  for c=1:Ci
    m_output{c} = zeros(size(m_input{c}));
  end
  % add the codons from each sampled block to the mask (some of them more than once!)
  for j=1:length(iblocks)
    iC   = block_chr(iblocks(j));
    pos_block_start = 1 + cfg.bootstrap.block_size * (block_idx(iblocks(j))-1);
    iidx = find( pos_input{iC} >= pos_block_start & pos_input{iC} < pos_block_start+cfg.bootstrap.block_size );
    m_output{iC}(iidx) = m_output{iC}(iidx) + m_input{iC}(iidx);
  end
  
  % write to file
  for c=1:Ci
    f = fopen( sprintf('%s_%03d.txt', files_bootstrap.output_mask_pref{c}, i), 'wt' );
    for k=1:length(pos_input{c})
    fprintf(f, '%d\t%d\n', pos_input{c}(k), m_output{c}(k) );
    end
    f = fclose(f);
  end
end



% % mode 2: split the dataset to two halfs
% TBD


end

