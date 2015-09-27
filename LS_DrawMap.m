
function [DivRedPred, stats] = LS_DrawMap( outmap_pref, params, GEs, positions, cfg_inf )


if isstr(cfg_inf)
  cfg_inf = file2struct(cfg_inf);
end

if isstr(GEs)
  files_buildGE_file = GEs;
  GEs    = LS_PrecalcGridElements(  files_buildGE_file, cfg_inf.GEs );
end


if isempty(positions) | (isstr(positions) & strcmp(positions,''))
  for c=1:length(cfg_inf.output.chromosomes)
    if isfield(GEs,'BSbase') & ~isempty(GEs.BSbase)
      chr_len = GEs.BSbase{c,1}.cfg.chr_len;
    else
      chr_len = GEs.SWbase{c,1}.cfg.chr_len;
    end
    positions{c} = unique([[1:cfg_inf.output.spatial_resolution:chr_len] chr_len]);
  end
elseif isstr(positions)
  % load for each chromosome
  spositions = file2struct(positions);
  positions{c} = load( spositions.positions{c} );
end

% a bit patchy...
tau       = params(3);
%

if isfield(cfg_inf, 'theta0')
  EgMutDiv  = cfg_inf.theta0*tau;
else
  EgMutDiv  = 0;
end

C        = length(cfg_inf.output.chromosomes);
nSWannos = size( GEs.SWbase, 2 );
nBSannos = size( GEs.BSbase, 2 );

for c=cfg_inf.output.chromosomes
    
  [~, ~, stats.SWparams{c}, stats.BSparams{c}, DivRedPred{c}] = composePredictedDiversity2( GEs.SWbase(c,:), GEs.gFocGrid{c}, GEs.BSbase(c,:), params, cfg_inf.inf, positions{c}, EgMutDiv );
  
%   [~, ~, stats.SWparams{c}, stats.BSparams{c}, DivRedPred{c}] = composePredictedDiversity2( GEs.SWbase(c:nSWannos:end), GEs.gFocGrid{c}, GEs.BSbase(c:nBSannos:end), params, cfg_inf.inf, positions{c}, EgMutDiv );
  
  % write to file
  SaveLSMap( [outmap_pref cfg_inf.chr_id{c}], DivRedPred{c}, cfg_inf.output.LSmap_res );
  
end


end

