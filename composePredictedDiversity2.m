
function [gSWj, gBSj, SWparams, BSparams, DivRedPred] = composePredictedDiversity2( SWbase, SWgpos, BSbase, params, config, pos, EgMutDiv ) %, remote2closeDiv
% function [DivRedPred, sumeps] = composePredictedDiversity2( preLH, params, config, gpos )

% calculating the predicted reduction in diversity due to linked selection
% (selective sweeps and background selection)
% redirects from calcPreLHstats initially

global MLParamsStruct;

pos = double(pos);

% prepare SW and BS grid elements for all spatial positions
% fill the thinned grid in to reincorporate all positions for this stage of
% the inference


gSWj = [];
if ~isempty(SWbase{1}) % as long as at least 1 anno exists SWbase{1} will not be empty
  for a=1:length(SWbase) % length of annos (2 right now)
    for k=1:length(SWbase{a}.gSWj) % length of selection coef vector
      idx = ~isnan(SWbase{a}.gSWj{k}); % idx all non-nan vals (all)
      gSWj{a}(k,:) = interp1(double(SWgpos.pos{k}(idx)), SWbase{a}.gSWj{k}(idx), pos(:)'); % interpolate out the selection on grid points to surrounding sites to use all neutral sites
      if config.use_fake_subs % is NOT used on human run...
        gSWj{a}(k,:) = gSWj{a}(k,:) + interp1(double(SWgpos.pos{k}(idx)), SWbase{a}.gSWj_fake{k}(idx), pos(:)');
      end
      gSWj{a}(k,isnan(gSWj{a}(k,:))) = 0;  % sites that are NaN after the interp1, turn to 0.
    end
  end
else
  gSWj{1} = ones(size(pos))';
end

% this process for bs is the same as that for sw but now instead of coal
% rates each position has some <= 0 selection coef 
gBSj = [];
if ~isempty(BSbase{1})
  for a=1:length(BSbase)
    for k=1:length(BSbase{a}.Bj)
      cX = cumsum([0 BSbase{a}.Lj{k}(1:end-1)]); % collect the ranges of BS effects and order them cumulatively
      cX = [0.75+cX; BSbase{a}.Lj{k}+0.25+cX]; % turn these into segments for interpolation
      cB = double([BSbase{a}.Bj{k} BSbase{a}.Bj{k}])'; % values for background selection
      gBSj{a}(k,:) = interp1(double(cX(:)), double(cB(:)), double(pos));
    end
  end
else
  gBSj{1} = ones(size(pos))';
end


% collect (selection) parameters of the grid elements

SWparams = [];
if ~isempty(SWbase{1})
    
  for a=1:length(SWbase)
    SWparams{a}.s    = SWbase{a}.cfg.FE_grid; % s values
    SWparams{a}.Ne0  = SWbase{a}.cfg.Ne0; %SWbase{a}.params{k}.Ne0;
    SWparams{a}.remote2closeDiv = config.remote2closeDiv;
  end
  
else  % some defaults...
  SWparams{1}.s = 1;
  SWparams{1}.Ne0 = 1;
  SWparams{1}.remote2closeDiv = 1;
end

BSparams = [];
if ~isempty(BSbase{1})
    
  for a=1:length(BSbase)
    BSparams{a}.u         = BSbase{a}.cfg.u_del; % del mut rate
    BSparams{a}.t         = BSbase{a}.cfg.FE_grid; % grid of t vals for neg selection
    BSparams{a}.nsites    = BSbase{a}.cfg.anno_len; % number of sites under selection (sites at diff annotations)
  end
  
else % some defaults...
  BSparams{1}.u         = 10^-10;
  BSparams{1}.t         = 1;
  BSparams{1}.nsites    = 1;

end



% construct BS, SW and combined linked selection maps from the grid elements and parameters

if ~isempty(params) % this WILL be empty on the initial precalc that brings us here in the beginning of the LH maximization, it is just an empty cell []

  % re-write this to work with memmap versions of the two selection structs
  DivRedPred = composeLSMapFromElements( 0, gSWj, gBSj, params, SWparams, BSparams, EgMutDiv, config );
  
  DivRedPred.pos = pos;

else
  
  DivRedPred = [];  % if there are no parameters, there is no reduction prediction to make
  
end


end
