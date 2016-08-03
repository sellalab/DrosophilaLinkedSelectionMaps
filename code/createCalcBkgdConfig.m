
function errval = createCalcBkgdConfig( conf_file, config )

errval = 0;

% config.Et = 10^-4;
% config.u_del = 10^-7;
% config.t_dist_type = 'POINT';
% config.t_dist_gamma_shape = 0.3;
% 
% config.chr_features = 'E:\Downloads\human\chromosome_features.txt';
% config.rec_table = 'E:\Downloads\human\genetic_maps\Myers\genetic_map_chr1_b36.txt';
% config.rec_table_scale = 10^-8;
% config.cons_table = 'E:\Downloads\human\data_McVicker\gcons_ex\features\chr1.coords';
% 
% config.output_dir = 'E:\code\vc\bkgd\debug\test\';
% config.output_token = '_test';
% config.B_res = 10;


f = fopen(conf_file, 'wt');

fprintf(f, 'CHROMOSOME_NAME=%s\n',   config.chr_id); %E:\Downloads\human\chromosome_features.txt
fprintf(f, 'CHROMOSOME_LENGTH=%d\n', config.chr_len); %E:\Downloads\human\chromosome_features.txt

fprintf(f, 'CHROMOSOME_FEATURES=%s\n', config.chr_features); %E:\Downloads\human\chromosome_features.txt
fprintf(f, 'OUTPUT_DIR=%s\n', config.output_dir); %E:\code\vc\bkgd\debug\test\
fprintf(f, 'OUTPUT_TOKEN=%s\n', config.output_token); %_t0.0001_lrB
fprintf(f, 'RECOMB_RATE_TABLE=%s\n', config.rec_table); %E:\Downloads\human\genetic_maps\Myers\genetic_map_chr1_b36.txt
fprintf(f, 'RECOMB_RATE_SCALE=%g\n', config.rec_table_scale); %0.00000001
fprintf(f, 'CONS_TABLE=%s\n', config.cons_table); %E:\Downloads\human\data_McVicker\gcons_ex\features\chr1.coords
fprintf(f, 'PARAM_T=%g\n', config.Et); %0.0001
fprintf(f, 'PARAM_U=%g\n', config.u_del); %1e-7
fprintf(f, 'PARAM_T_DIST_TYPE=%s\n', config.t_dist_type); %POINT
if isfield(config, 't_dist_gamma_shape')
    fprintf(f, 'PARAM_GAMMA_SHAPE=%f\n', config.t_dist_gamma_shape); %0.3
end
fprintf(f, 'BKGD_SCALE=%f\n', config.B_res);% 10.0

fprintf(f, 'USE_SUM_APPROXIMATION=0\n');
fprintf(f, 'PARAM_T_DIST_TRUNC_LOW=0.00001\n');
fprintf(f, 'PARAM_T_DIST_TRUNC_HIGH=1.0\n');
fprintf(f, 'BKGD_INTERP_MAX_DIST=0.0001\n');
fprintf(f, 'BKGD_CHANGE_THRESH=0.02\n');
fprintf(f, 'BKGD_OFFSET=0.01\n');
fprintf(f, 'BKGD_PARAM_MAX_SUM_THRESH=0.001\n');

f = fclose(f);

end

