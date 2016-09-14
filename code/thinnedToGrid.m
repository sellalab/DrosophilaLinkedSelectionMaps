
function [thinned_pos, thinned_idx, needed] = thinnedToGrid( pos, grid )

[pos, idx] = unique(pos);

delta   = grid(1);
pos_min = grid(2);
pos_max = grid(3);


% needed = zeros(size(pos));
% 
% for curpos = pos_min:delta:pos_max
%   needed(find(pos>=curpos,1,'first')) = 1;
%   needed(find(pos<=curpos,1,'last' )) = 1;
% end
% thinned_pos = pos(find(needed));

% HISTC(X,EDGES)
% edit 12/27 create a nan filter to make sure we dump positions outside the
% scope of the gmap:
% NOTE: this maybe not be sufficient downstream when the positions come
% back, instead, alter the data or the gmaps themselves so that all neutral
% sites are within the gmap
% nanfilter = isfinite(pos);
% needed = histc([pos_min:delta:pos_max], [0; pos(nanfilter) ;pos_max + 1]);
needed = histc([pos_min:delta:pos_max], [0; pos ;pos_max + 1]);
thinned_i1 = find(needed(1:end-2)>0);
thinned_i2 = [thinned_i1(2:end)-1  length(pos)];

thinned_pos = pos([thinned_i1 thinned_i2]);
thinned_idx = idx([thinned_i1 thinned_i2]);

[thinned_pos, idx0] = unique(thinned_pos);
thinned_idx = thinned_idx(idx0);


bigo = 7;
