
function te = SaveSWBase( output_pref, gFocGrid, SWbase )

% write SW levels for a single chromosome for a single annotation

for t=1:length(gFocGrid.pos)
  
  te.outputfiles{t} = sprintf( '%s_t%.6f.sw', output_pref, SWbase.CalcSW.S(t) ); %SWbase.params{t}.S
  f = fopen( te.outputfiles{t}, 'wt' );
  
  % write header with parameters
  
  fprintf( f, '#%s %d\n',                     SWbase.CalcSW.chr_id, SWbase.CalcSW.chr_len );
  fprintf( f, '#CHROMOSOME_NAME=%s\n',        SWbase.CalcSW.chr_id );
  fprintf( f, '#CHROMOSOME_LENGTH=%d\n',      SWbase.CalcSW.chr_len );
 %fprintf( f, '#CHROMOSOME_FEATURES=%s\n',    SWbase.chr_features );
  fprintf( f, '#OUTPUT_FILE=%s\n',            te.outputfiles{t} );
  fprintf( f, '#OUTPUT_TOKEN=%s\n',           SWbase.CalcSW.output_token );
  fprintf( f, '#RECOMB_RATE_TABLE=%s\n',      SWbase.CalcSW.rec_table ); 
  fprintf( f, '#RECOMB_RATE_SCALE=%e\n',      SWbase.CalcSW.rec_table_scale ); %10^-8
  fprintf( f, '#FOCALS_TABLE=%s\n',           SWbase.CalcSW.focals_table ); 
  fprintf( f, '#PARAM_NE0=%e\n',              SWbase.CalcSW.Ne0 ); %SWbase.params{t}.Ne0
  fprintf( f, '#PARAM_S=%e\n',                SWbase.CalcSW.S(t) ); %SWbase.params{t}.S
  fprintf( f, '#PARAM_S_DIST_TYPE=%s\n',      'POINT' );
  fprintf( f, '#TTF_APPROXIMATION=%s\n',      SWbase.CalcSW.trap_aprx ); %'diffusion'
%   fprintf( f, '#USE_FAKE_SUBSTITUTIONS=%d\n', SWbase.CalcSW.use_fake_subs ); 
  fprintf( f, '#StopSum=%d\n',                SWbase.CalcSW.StopSum ); 
  fprintf( f, '#INTERP_METHOD=%s\n',          SWbase.CalcSW.InterpMethod ); 
  fprintf( f, '#MAX_DIST_SCALED=%d\n',        SWbase.CalcSW.gMaxDistScaled ); 
  fprintf( f, '#MAX_DIST=%e\n',               SWbase.CalcSW.gMaxDist(min(t, length(SWbase.CalcSW.gMaxDist))) ); 
  
  
  % write data - physical position and CSW coef values (genetic position is
  % relevant but redundant - can be calculated from the genetic map files
  % and save lots of space)
  if ~isfield(SWbase, 'gSWj_fake')
    SWbase.gSWj_fake{t} = zeros(size(SWbase.gSWj));
  end
  
  for i=1:length(gFocGrid.pos{t})
    fprintf( f, '%d\t%.3e\n', gFocGrid.pos{t}(i), SWbase.gSWj{t}(i));%, \t%.3e SWbase.gSWj_fake{t}(i) DM edit
  end
  
  f = fclose( f );
  
end

end

