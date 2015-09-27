
function [outpos, outrat] = applyGenmap2pos(genmap, inpos, invert)

if nargin < 3
  invert = 0;
end

if ~invert
%     'condition 1' THIS IS THE CONDITION THAT RUNS IN THE DEFAULT CONFIG
    % the point of this function is to interpolate the map distance between SW sites based on the genetic map we have
  outpos  = interp1( genmap.pos, genmap.R, double(inpos) ) / 100; % DM edit convert to Morgans...
  % c(i) = ( R(i)-R(i-1) ) / ( x(i)-x(i-1) ), so actually it corresponds to (x(i-1)+x(i))/2
  if nargout>1
%      'condition 2'
    outrat  = interp1( [genmap.pos(1); 0.5*(genmap.pos(1:end-1)+genmap.pos(2:end))], genmap.c, double(inpos) );
%     outrat  = interp1( genmap.pos, genmap.c, double(inpos) );
  end
else
%    'condition 3'
  outpos  = interp1_weakmonotonic(genmap.R/100, genmap.pos, double(inpos));
%   i1 = find(~isnan(outpos), 1, 'first');
%   i2 = find(~isnan(outpos), 1, 'last');
%   outpos(1:i1-1)   = outpos(i1);
%   outpos(i2+1:end) = outpos(i2);
  outrat  = [];
end

end

