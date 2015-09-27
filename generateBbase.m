
function [config, Bbase] = generateBbase(config, t, skip_generate)

if nargin < 3
  skip_generate = 0;
  'skipping generation of new Bmaps'
end

%from precalc function: generateBbase(CalcBS, CalcBS.FE_grid, CalcBS.skip_generate_maps);


% % params & configuration
% 
% config.L = 1000000; %chr len in bp
% config.Et = 10^-4;
% config.u_del = 10^-8;
% config.t_dist_type = 'POINT';
% config.N4u = 0.001;
% config.T2u = 0.01;
% 
% config.chr_features = 'E:\projects\human\simulated\chromosome_features.txt';
% config.rec_table = 'E:\projects\human\simulated\genetic_map_simulated.txt';
% config.rec_table_scale = 10^-8;
% config.cons_table = 'E:\projects\human\simulated\cons.coords';
% config.mutvar_map = 'E:\projects\human\simulated\';
% 
% config.output_dir = 'E:\projects\human\simulated\B\';
% config.output_token = '_test';
% config.B_res = 1000;


if isempty(config)
    % not empty on the example run
    'making defaultFakeBSdataConfig'
    config = defaultFakeBSdataConfig();
end

if ~isfield(config, 'output_pref')
    % why set this to empty string?
  config.output_pref = '';
end

%% load features
% required features are read within 'createCalcBkgdConfig'

% create B maps

% exefile = 'E:\projects\sweeps\NeutFoc\BS_util\BSmap.exe'; %'E:\projects\human\simulated\calc_bkgd.exe';
% exefile = 'E:\code\vc\bkgd\debug\BSmap.exe'; %'E:\projects\human\simulated\calc_bkgd.exe';
% exefile = 'E:\projects\sweeps\NeutFoc\BS_util\calc_bkgd.exe'; %'E:\projects\human\simulated\calc_bkgd.exe';

% empty string for pref is just an optional filename suffix
output_dir = [config.output_dir config.output_pref];

config.t_dist_type = 'POINT';
cfgfilename  = '';


unique_token = sprintf('%s_t%f', config.name, t(1));
for i=2:length(t)
  unique_token = sprintf('%s_t%f', unique_token, t(i));
end
if isfield(config, 'unique_id')
  unique_token = sprintf('%s_%s', unique_token, config.unique_id);
end
config.unique_id = '';


for k=1:length(t)
%   'entering createCalcBkgdConfig loop'
    config.Et = t(k);
    config.output_token = sprintf('%s_t%f', config.name, config.Et);
    config.basefiles{k} = sprintf('%s.bkgd', config.output_token); %config.chr_id, 
    config.grid_files{k}= [output_dir config.basefiles{k} '.bkgd']; % something wrong here? I dont see any .bkgd.bkgd files around, maybe an error?

    curdiri = cd();
% %     while isempty(cfgfilename) | exist(cfgfilename)==2
%       randizi = randi(10^8, 1);
% %       cfgfilename  = sprintf('%s/simuBconf_%d.txt', curdiri, randizi); 
% %     end
%     cfgfilename2 = sprintf('%s/simuBconf_%d.txt', curdiri, randizi); 
    cfgfilename = sprintf('%s/simuBconf_%s.txt', curdiri, unique_token); 

    createCalcBkgdConfig( cfgfilename, config );
    comman = sprintf( '!%s %s %s', config.exe_file, cfgfilename, config.output_token ); 
    if skip_generate==0
      T = evalc(comman) % DM edit: I want to see the printout so un-supress here
    end
    
    delete(cfgfilename); %DM edit- keep this to check out 

end


% read B map
if nargout>1
  Bbase = loadBbase(config.basefiles, output_dir);
end


end

