
function [vEpsilon, vS, vS_weights] = trapProbability(vR, Ne, vS, vS_weights, sweep_approx)

% this function reproduces Durrett's book tables, with input N/2

eulero = -psi(1);

L = length(vR);

if 0 %length(vR>5000)  % this is completely deactivated...
%   delta = 1/100;
%   % delta = 1/100;
%   towardsTail  = 2;
%   towardsFocal = 1;
%   
%   interp_eps = 1;
%   % subsample the radii vector, depending whether it is one-sided or two-sided
%   minEpsilon = max(10^-200, exp( -2*vR(end)./vS.*log(2*Ne) ));
%   C = -log(minEpsilon)./vR(end);
%   sub_vR = [vR(1) vR(end)];
%   for i=1:length(vS)
%     sub_vR = [sub_vR  -log(minEpsilon(i)+[0:2*delta:1].^towardsFocal.*(1-minEpsilon(i)))/C(i)  -log(minEpsilon(i)+[0:2*delta:1].^towardsTail.*(1-minEpsilon(i)))/C(i)];
%   end
%   sub_vR = sort(sub_vR);
else
    
  interp_eps = 0;  % don't change any of the inputs here...
  sub_vR = vR;
  
end
sub_vEpsilon = zeros(size(sub_vR));  % an empty vector for the results

for i=1:length(vS)
    
    
  if vS>0
      
    switch sweep_approx
      case 'deterministic'
        tfix = timeToFixation( Ne, vS(i), 'deterministic' );
      
        
        
      case 'diffusion'  % use diffusion for the inference currently
        tfix = timeToFixation( Ne, vS(i), 'diffusion_approximation' );
%         tfix = timeToFixation( Ne, vS(i), 'Wakeley' );
%       case 'durrett'
%         vM = floor(4*Ne*vS); % M=2*Ne*s for haploids. Replace N->2N here and in the logistic model for diploids
%         % there must be some error in one of the approximations - the 4-factor used here should be 2 according to my derivation, but then the results do not fit the logistic sweep approximation
%         % Durrett, prop 3
%         % exp(-2*r/s*(log(M+1)+eulero-3/2))+2*r/s.*(1/2-1/(M+1)).*(M+1).^-(2*r/s) for haploids. Replace N->2N here and in M for diploids
%         coalprob = exp( -2*sub_vR/vS(i)*(log(vM(i)+1)+eulero-3/2))+2*sub_vR/vS(i).*(1/2-1/(vM(i)+1)).*(vM(i)+1).^-(2*sub_vR/vS(i));
    end
    
    % what is the point of this, it cancels out to 1?
    IntOneMinusFreq = 1/2; % ~(p1^2-p0^2)/2; % OLD BUGGED VALUE = 1;
    coalprob = exp( -2*sub_vR*tfix*IntOneMinusFreq );  % probability of coalescing for a given site based on its distance from substitution
    
    sub_vEpsilon = sub_vEpsilon + vS_weights(i)*coalprob;
  end
end

if interp_eps  % FALSE for now, this section is skipped...
  vEpsilon = interp1q(sub_vR',sub_vEpsilon',vR); % interp1(sub_vR,sub_vEpsilon,vR, 'spline');
  if isnan(vEpsilon(end))
    vEpsilon(end) = vEpsilon(end-1);
  end
  if(sum(isnan(vEpsilon)))
    bigo = 7;
  end
else
  vEpsilon = sub_vEpsilon;
end

end

