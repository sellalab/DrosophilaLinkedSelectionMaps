function config = defaultFakeBSdataConfig()

% params & configuration

config.L = 1000000; %chr len in bp
config.chr_id = 'chr1';

% this is the actual mean selection strength that will be used by calc_bkgd.exe
config.Et = 10^-4;
% these vectors are used by higher functions that may want to create B maps
% that are combinations of basic B-atoms. Each component is created by the
% higher function by substitutiong for each i Et=vEt(i) thus producing a
% base, from which the final Bmap is composed
config.vEt = 10^-4; %[10^-2 10^-3 10^-5]
config.vWt = 1;     %[0.3   0.1   0.6  ]

config.u_del = 10^-8;
config.t_dist_type = 'POINT'; % GAMMA / POINT / EXPONENTIAL
config.N4u = 0.001*10;
config.T2u = 0.01*10;

config.chr_features = 'E:\projects\human\simulated\chromosome_features.txt';
config.rec_table = 'E:\projects\human\simulated\genetic_map_simulated.txt';
config.rec_table_scale = 10^-8;
config.cons_table = 'E:\projects\human\simulated\cons.coords';
config.mutvar_map = 'E:\projects\human\simulated\';

config.output_dir = 'E:\projects\human\simulated\B\';
config.output_token = '_test';
config.B_res = 1000;

end

