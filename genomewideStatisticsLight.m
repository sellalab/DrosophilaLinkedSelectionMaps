
function gwStats = genomewideStatisticsLight( fdata )


if nargin < 2
  poolNcorrStats = 0;
end


C = length(fdata);

n_poly_sites = 0; 
N_sites = 0; 
sum_pq_product = 0; 
comparisons = 0;

gMutDiv = [];
gwStats.mapL = 0;

for c=1:C
  
  % this is the only part of the script that is used now. the labels are
  % kind of misleading, make better labels in python...
  % calculation of pi = p * q / (N choose 2)
  
  n_poly_sites = n_poly_sites + length(fdata{c}.Poly.idxS); % number of polymorphic sites
  sum_pq_product = sum_pq_product+ sum(prod(fdata{c}.Poly.specS,2)); % p*q for calculating het (counts of allele A and B)
  N_sites = N_sites + length(fdata{c}.Poly.gpos); % total alleles
  comparisons= comparisons+ sum(fdata{c}.Poly.sample.*(fdata{c}.Poly.sample-1)/2); % heterozygosity (pi)
  
  
  pos = double(fdata{c}.Poly.pos);
  dSidx = find( ~isnan(fdata{c}.MutProx.dS) | fdata{c}.MutProx.dS < 0 ); 
  cgMutDiv   = interp1(double(fdata{c}.MutProx.pos(dSidx)), fdata{c}.MutProx.dS(dSidx), pos); % interpolates mut rate (only two points needed for constant)
  cgMutDiv(pos < fdata{c}.MutProx.pos(dSidx(1  ))) = fdata{c}.MutProx.dS(dSidx(1  ));
  cgMutDiv(pos > fdata{c}.MutProx.pos(dSidx(end))) = fdata{c}.MutProx.dS(dSidx(end));
  gMutDiv = [gMutDiv; cgMutDiv];
  
  gwStats.mapL = gwStats.mapL + max(fdata{c}.Poly.gpos) - min(fdata{c}.Poly.gpos); % get the full genetic map length in Morgans
  gwStats.chr_len(c)       = fdata{c}.chr_len;
  
end

gwStats.Poly            = (n_poly_sites/N_sites); % this is the statistic for average polymorphism
gwStats.Het             = (sum_pq_product/comparisons); % this is the statistic for average diversity
gwStats.MutProx         = nanmean(gMutDiv); % this is the statistic for the average mutation rate
  

end

