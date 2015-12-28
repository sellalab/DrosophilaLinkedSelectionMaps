%% loading a base of atom-based B maps

function Bbase = loadBbase(maps, mapsdir)

if nargin==1 | isempty(mapsdir)
    mapsdir = '';
% else
%     mapsdir = [mapsdir '\'];
end

% mapsdir = 'E:\code\vc\bkgd\debug\test';
% Bbase.t = [.99 .1 .01 .001 .0001];
% Bbase.res = [1000 1000 1000 1000];

for k=1:length(maps)
	% maps{k}  -- debug printout (12/18) to see that the map names are formatted correctly
    [BL{k}, cfg{k}] = loadBmap( [mapsdir maps{k}] );
%    [BL{k}, params] = loadBmap( [mapsdir '\' sprintf('chr1_ex_t%g.bkgd', Bbase.t(k))] );
end



for k=1:length(BL)
  rpos = cumsum(BL{k}(:,2));
  Bbase.Lj{k} = [rpos(1); diff(rpos)]';
  Bbase.Bj{k} = log(double(BL{k}(:,1))/cfg{k}.res); % WHY?
end
Bbase.cfg = cfg;

