

function xx = epsilonFromBoundary( x, rangeL, rangeH, rel_epsilon )

xx = x;


ii = find( x==rangeL & rangeL<rangeH );

xx(ii) = rangeL(ii)*(1-rel_epsilon) + rangeH(ii)*rel_epsilon;

ii = find( x==rangeH & rangeL<rangeH );

xx(ii) = rangeH(ii)*(1-rel_epsilon) + rangeL(ii)*rel_epsilon;


end

