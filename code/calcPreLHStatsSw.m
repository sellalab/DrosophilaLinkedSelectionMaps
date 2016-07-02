
function preLHstats = calcPreLHStatsSw( SWbase, SWgpos, BSbase, data, config, weights, ccounts, pos_external, gpos_external ) %remote2closeDiv, 


% min_paml_codons = 100;
% scale_dS = 1; %0.4;


if ~exist('gpos_external')
  preLHstats.gpos_ext = [];
else
  preLHstats.gpos_ext = gpos_external;
end
if ~exist('pos_external')
  preLHstats.pos_ext = [];
else
  preLHstats.pos_ext = pos_external;
end
if ~exist('ccounts')
  ccounts = [];
end
if ~exist('weights')
  weights = ones(size(data.Poly.pos)); % every position given a weight of 1 to start with
end
if length(weights) ~= length(data.Poly.pos)
  error('predef_idx length different from Poly.pos length')
end

%% apply masks if masks are used


idx1p = data.Poly.idxS; % this is the index of the poly sites in the vector of poly sites (not real pos)
preLHstats.idxc = find(weights>0); % this is the size of the entire polymorphism set since no masking
[te, preLHstats.idxc1p, iidx1p] = intersect(preLHstats.idxc, idx1p); % here segregating sites that pass masks are taken out (all for human)
preLHstats.idxc1p = preLHstats.idxc1p'; % index 1p means there is a polymorphism (polymorphism==TRUE), also take the inversion for format reasons
preLHstats.gpos = data.Poly.gpos(preLHstats.idxc); % get genetic position on the indices we have kept from above filters (everything is kept)
preLHstats.pos  = data.Poly.pos( preLHstats.idxc); % get physical positions as well


%% calculate basic stats for the polymorphism dataset:


samples = [weights(preLHstats.idxc).*double(data.Poly.sample(preLHstats.idxc)) zeros([length(preLHstats.idxc) 1]) ]; % empty matrix of zeros to init
samples(preLHstats.idxc1p,:) = (weights(preLHstats.idxc(preLHstats.idxc1p))*[1 1]).*double(data.Poly.specS(iidx1p,:));

preLHstats.samplesHet = uint16(prod(samples, 2)); % pq for every site in the data
preLHstats.samplesHom = uint16(sum(samples,2).*(sum(samples,2)-1)/2 - double(preLHstats.samplesHet)); % N choose 2 - pq at every site
preLHstats.samplesHetSum = sum(preLHstats.samplesHet); % sum of the pq product
preLHstats.samplesHomSum = sum(preLHstats.samplesHom); % sum of the (N choose 2) - pq product
preLHstats.samplesSum    = preLHstats.samplesHomSum+preLHstats.samplesHetSum; % (N choose 2) * number of sites

dSidx = find( ~isnan(data.MutProx.dS) & data.MutProx.dS >= 0 ); %find(data.MutProx.paml_codons >= min_paml_codons);
preLHstats.gMutDiv   = interp1(double(data.MutProx.pos(dSidx)), data.MutProx.dS(dSidx), double(preLHstats.pos)); %data.MutProx.dS(dSidx)*scale_dS
preLHstats.gMutDiv(preLHstats.pos < data.MutProx.pos(dSidx(1  ))) = data.MutProx.dS(dSidx(1  )); %*scale_dS;
preLHstats.gMutDiv(preLHstats.pos > data.MutProx.pos(dSidx(end))) = data.MutProx.dS(dSidx(end)); %*scale_dS;

preLHstats.SumMutProx = sum(preLHstats.gMutDiv.*(double(preLHstats.samplesHet)+double(preLHstats.samplesHom)));
preLHstats.EgMutDiv = nanmean(preLHstats.gMutDiv);
preLHstats.EgHet    = preLHstats.samplesHetSum/(preLHstats.samplesSum); % average level of heterozygosity

%% compose predicted diversity


if isempty(preLHstats.pos_ext) | isempty(preLHstats.gpos_ext)
  cur_pos  = double(data.Poly.pos(preLHstats.idxc)); % the data just assembled will be used, we don't use external pos or gpos with the current settings
else
  cur_pos  = preLHstats.pos_ext;
end




[preLHstats.gSWj, preLHstats.gBSj, preLHstats.SWparams, preLHstats.BSparams] = composePredictedDiversity2(SWbase, SWgpos, BSbase, [], config, cur_pos, preLHstats.EgMutDiv);


bigo = 7;
