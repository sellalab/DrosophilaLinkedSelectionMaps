%% loading a base of atom-based B maps

function Bbase = loadBbase(maps, mapsdir)

if nargin==1 || isempty(mapsdir)
    mapsdir = '';
end

for k=1:length(maps)
    [BL{k}, cfg{k}] = loadBmap([mapsdir maps{k}]);
end



for k=1:length(BL)
  rpos = cumsum(BL{k}(:,2));
  Bbase.Lj{k} = [rpos(1); diff(rpos)]';
  Bbase.Bj{k} = log(double(BL{k}(:,1))/cfg{k}.res);
end
Bbase.cfg = cfg;

