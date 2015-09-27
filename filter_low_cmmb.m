function kpidx = filter_low_cmmb(genmap, pos)

% a function that scans the genetic map at each position and retrieves the
% recombination rate column

    % filter low recombination block:
    [~, ~, bins] = histcounts(pos, genmap.pos);
    filtered_bins = bins(bins > 0);
    cmmb = genmap.c(filtered_bins-1);
    kpidx = cmmb > 0;  % keep index
%     positions = pos(kpidx);