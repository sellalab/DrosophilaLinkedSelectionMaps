
function  SwCoef = SwCoef1point( S, Ne0, gFocSite, gSwSites, config )

    if config.StopSum == 2
        % complete...
        SwCoef = 0;
    else
        if config.StopSum == 1
            idx = find(abs(gSwSites-gFocSite)<config.gMaxDist);  % this is the setting that is used
        else
            idx = [1:length(gSwSites)];
        end
        % essentially returns a sum of the probabilities of fixation for
        % all sites close to putative sweep sites (within s Morgans of the
        % site)
        SwCoef = sum(trapProbability(abs(gSwSites(idx)-gFocSite), Ne0, S, 1, config.trap_aprx));
    end

end

