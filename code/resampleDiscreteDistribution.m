
function [w1, x1] = resampleDiscreteDistribution( x0, w0, x1, support_edges )

[x0, idx] = sort(x0);
w0        = w0(idx);

x1        = sort(x1);

% w0 = w0 / sum(w0);


% this is a grid of interpolation points for the continuous approximation of the cdf
x0i     = [support_edges(1) (x0(1:end-1)+x0(2:end))/2 support_edges(2)];
cdf0i   = [0 cumsum(w0)];

% edges of bins centered around x1 points
x1i     = [support_edges(1) (x1(1:end-1)+x1(2:end))/2 support_edges(2)];
cdf1i   = interp1( x0i, cdf0i, x1i, 'linear' );

w1      = diff( cdf1i );

if sum(w1)~=1
  bigo = 7;
end
