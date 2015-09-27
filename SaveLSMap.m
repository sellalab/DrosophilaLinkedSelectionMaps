
function errval = SaveLSMap( outpref, LSDivPreds, res )

errval = 0;

if nargin<3
  res = 100;
end


% pos, LS, BS, SW

  f = fopen( [outpref '.LS'], 'wt' );
  
  fprintf(f, '#RES=%d\n', res);
  fprintf(f, '#pos\tLS\tBS\t1/(1+CS)\n');
  
  for i=1:length(LSDivPreds.pos)
    fprintf(f, '%d\t%d\t%d\t%d\n', LSDivPreds.pos(i), round(res*LSDivPreds.Red(i)), round(res*LSDivPreds.BS(i)), round(res*1/(1+LSDivPreds.cSW(i))));
  end
    
  f = fclose(f);

end

