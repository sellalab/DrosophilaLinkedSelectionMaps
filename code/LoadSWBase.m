
function [SWbase, gFocGrid] = LoadSWBase( input_pref, cfg )


coefs = cfg.FE_grid;

for t=1:length(coefs)

  te.inputfiles{t} = sprintf( '%s_t%.6f.sw', input_pref, coefs(t) );

  [SWbase.gSWj{t}, SWbase.gSWj_fake{t}, gFocGrid.pos{t}, cfg] = loadSWmap( te.inputfiles{t}, cfg );

  % we assume that all configuraion params besides the following fields are the same for all grid elements, therefore
  % cummulating these to assign their list of values at the end of element loading
  grid_files{t} = cfg.output_file;
  FE_grid(t)    = cfg.S;
  gMaxDist(t)   = cfg.gMaxDist;

end

cfg.grid_files = grid_files;
cfg.FE_grid    = FE_grid;
if length(unique(gMaxDist))>1
  cfg.gMaxDist = gMaxDist;
end

SWbase.cfg = cfg;
