
function vSw = SwCoef( vS, Ne0, gFocSites, gSwSites, cfg ) %md_correction_used

% arguments in LS_PrecalcGridElements: 
% SWbase{c,b} = SwCoef( CalcSW.FE_grid, CalcSW.Ne0, gFocGrid{c}, {Annots.SW{c}{b}.focals.gpos(Annots.SW{c}{b}.focals.isfake==0)}, CalcSW );
% vS - vector of deleterious fitness effect sizes
% Ne0 - effective population size, prior to the effects of sweeps
% gFocSites - positions of neutral sites where the SW coefficient is calculated, in genetic distance.
%   This may be only a grid of focal sites, and exact values for other sites can be interpolated at higher levels (beyond this function)
% gSwSites - positions of potential sweeps, in genetic distance (usually amino-acid substitution sites)


% DURRETT1 = 1;  DURRETT2 = 2;


if ~isfield(cfg, 'StopSum '),         cfg.StopSum       = 1;              ,end
if ~isfield(cfg, 'InterpMethod '),    cfg.InterpMethod  = 'linear';       ,end
if ~isfield(cfg, 'trap_aprx '),       cfg.trap_aprx     = 'diffusion';    ,end  
if ~isfield(cfg, 'gMaxDistScaled'),   cfg.gMaxDistScaled= 1;              ,end  
if ~isfield(cfg, 'gMaxDist'),         cfg.gMaxDist      = 1;              ,end  

gMaxDist = cfg.gMaxDist;

cfg.Ne0 = Ne0;
cfg.S   = vS;

for k=1:length(vS)
  
  vSw.gSWj{k} = zeros([length(gFocSites.gpos{k}) 1]);
  for i=1:length(gFocSites.gpos{k})
    if cfg.gMaxDistScaled
      cfg.gMaxDist = vS(k)*gMaxDist;  %DM currently used setting
    else
      cfg.gMaxDist = gMaxDist(k);
    end
    vSw.gSWj{k}(i) = SwCoef1point( vS(k), Ne0, gFocSites.gpos{k}(i), gSwSites{min(k,length(gSwSites))}, cfg );
  end
  if sum(isnan(vSw.gSWj{k}))>0
    eidx0              = find( isnan(vSw.gSWj{k}));
    eidx1              = find(~isnan(vSw.gSWj{k}));
    vSw.gSWj{k}(eidx0) = interp1(gFocSites.gpos{k}(eidx1), vSw.gSWj{k}(eidx1), gFocSites.gpos{k}(eidx0));
  end
  
  % a bit redundant - all information is in 'cfg'
  vSw.params{k}.Ne0 = Ne0;
  vSw.params{k}.S   = vS(k);
  
end


% vSw.md_correction = md_correction_used;
cfg.gMaxDist = gMaxDist;
vSw.CalcSW = cfg;

end
