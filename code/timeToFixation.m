
function ttf = timeToFixation( Ne, s, method, p0, p1 )

% approximating the time to fixation of a favorable allele with heteryzygous selective
% advantage s in a population of N diploids from a single copy (p0=1/2N)

if ~exist('p0')
  p0 = 1/(2*Ne);
end
if ~exist('p1')
  p1 = 1 - 1/(2*Ne);
end

switch method
  case 'diffusion_exact' % Ewens 1979 5.51+5.52 (note that formula 5.51 contains an error in this edition)
    ttf = 2*Ne*integral(@(x)t_condfix_p1(x, p0, 4*Ne*s), p0, p1); % 2s is the selective advantage of the homozygous, to be consistent with Ewens79 definitions
    
  case 'diffusion_approximation' % Hermisson and Pennings 2005 (appendix)
    eulero = -psi(1);
    ttf = 2*(log(4*Ne*s)+eulero-1/(4*Ne*s))/s;
    
  case 'deterministic' % logistic growth, starting from p=1/2N
    ttf = 2*log(2*Ne)/s;
    
end

end

%%

% function tt = t_condfix_0p(x,p,a)
% % 0<=x<=p
% % tt = 2*exp(-a*x).*exp(-a).*(1-exp(a*(1-p))).*(exp(a*x)-1).^2 ./ (a*x.*(1-x).*(1-exp(-a)).*(exp(-a*p)-1));
%
% logexpA_1 = log(exp(A)-1);
% if isnan(logexpA_1) | isinf(logexpA_1)
%   logexpA_1 = A;
% end
%
%
% log1_expA1_p          = log(exp(A*x)-1);
% logexpAx
% log1_exp_A
% logexp_Ap_1
%
% tt = exp( log(2) - a*x - a + log(1-exp(a*(1-p))) + 2*log(exp(a*x)-1) - log(a) - log(x.*(1-x)) - log(1-exp(-a)) - log(exp(-a*p)-1) );
% % tt = 2*exp(-a*x).*exp(-a).*(1-exp(a*(1-p))).*(exp(a*x)-1).^2 ./ (a*x.*(1-x).*(1-exp(-a)).*(exp(-a*p)-1));
%
% end

%%

function tt = t_condfix_p1(x,p,A)
% p<=x<=1
% tt = 2*(exp(A*x)-1).*(exp(A*(1-x))-1)                     ./ (A.*x.*(1-x).*(exp(A)-1));

logexpA_1 = log(exp(A)-1);
if isnan(logexpA_1) | isinf(logexpA_1)
  logexpA_1 = A;
end

logexpAx_1          = log(exp(A*x)-1);
logexpAx_1(A.*x>100)   = A.*x(A.*x>100);
logexpA1_x_1        = log(exp(A*(1-x))-1);
logexpA1_x_1(A*(1-x)>100) = A*(1-x(A*(1-x)>100));

tt = exp( log(2) + logexpAx_1 + logexpA1_x_1 - log(A) - log(x.*(1-x)) - logexpA_1);
% tt = 2*(exp(A*x)-1).*(exp(A*(1-x))-1)                     ./ (A.*x.*(1-x).*(exp(A)-1));

end
