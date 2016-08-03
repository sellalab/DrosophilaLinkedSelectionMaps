
function [R2, Dev2_obs_pred, Dev2_obs_Eobs] = Rsquare( obs, pred, weights, equalizemeans )

if nargin < 3 | isempty(weights)
  weights = ones(size(obs));
end

if nargin < 4
  equalizemeans = 0;
end


if equalizemeans
  pred = pred/sum(weights.*pred)*sum(weights.*obs);
end

Eobs = sum(weights.*obs)/sum(weights); %mean(obs);

Dev2_obs_pred = sum(weights.*(obs-pred).^2)/sum(weights);
Dev2_obs_Eobs = sum(weights.*(obs-Eobs).^2)/sum(weights);

R2 = 1 - Dev2_obs_pred/Dev2_obs_Eobs;

R2 = max(R2,0);

end

