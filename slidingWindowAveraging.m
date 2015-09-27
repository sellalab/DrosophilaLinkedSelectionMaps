
function sy = slidingWindowAveraging( x, y, window_size, window_jump, minimal_samples_per_window )

if ~exist('minimal_samples_per_window')
  minimal_samples_per_window = 1;
end

y = double(y);

if isempty(x)
  x = [0:length(y)-1];
else
  [x,idx] = sort(double(x));
  x = x-x(1);
  y = y(idx);
end


% wins = (x(end)-x(1))/win_jump + 1;

cwindow_size = window_size/window_jump;
cx = double(x)/window_jump+0.5-1;
% cx = floor((x+win_jump/2)/win_jump);
% rcx = round(cx);
% ucx = unique(rcx);



% patching an extra bin for the last incomplete window, to avoid possible
% extrapolation at the end of the grid
if floor(cx(end))==cx(end)
  wx =  [0:floor(cx(end))];
else
  wx = [[0:floor(cx(end))] cx(end)];
end

wy = NaN*ones(size(wx));

for i=1:length(wx)
  
  cidx = find( cx>=wx(i)-0.5*cwindow_size & cx<wx(i)+0.5*cwindow_size );
  
  window_samples = length(cidx);
  
  if window_samples > minimal_samples_per_window
    wy(i) = sum(y(cidx)) / window_samples;
  end
  
end



if sum( isnan(wy([1 end])) ) %sum( isnan(cy) )
  error('not enough samples in boundary windows'); % error('not enough samples in some windows');
end
  
sy = interp1( wx(~isnan(wy)), wy(~isnan(wy)), cx, 'pchip' );


