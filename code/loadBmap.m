
function [BL, cfg] = loadBmap( file )

% extract parameters from the commented lines

% cfg.res = 1000;
% file

f = fopen(file,'rt');

if f==-1
  error('BS-map file not found');
end

line = fgetl(f);

[cfg.chr_id, cfg.chr_len] = sscanf(line, '%s %d');

while (~isempty(line) & (line~=-1) & (line(1)=='#'))
    [x, y] = strtok(line(2:end), '=');
    switch x
        case 'CHROMOSOME_FEATURES'
            cfg.chr_features = y(2:end);
        case 'CHROMOSOME_NAME'
            cfg.chr_id = y(2:end);
        case 'CHROMOSOME_LENGTH'
            cfg.chr_len = str2num(y(2:end));
        case 'OUTPUT_DIR'
            cfg.output_dir = y(2:end);
        case 'OUTPUT_TOKEN'
            cfg.output_token = y(2:end);
        case 'RECOMB_RATE_TABLE'
            cfg.rec_table = y(2:end);
        case 'RECOMB_RATE_SCALE'
            cfg.rec_table_scale = str2num(y(2:end));
        case 'CONS_TABLE'
            cfg.cons_table = y(2:end);
        case 'PARAM_T'
            cfg.t = str2num(y(2:end));
        case 'PARAM_U'
            cfg.u = str2num(y(2:end));
        case 'PARAM_T_DIST_TYPE'
            cfg.t_dist_type = y(2:end);
        case 'PARAM_GAMMA_SHAPE'
            cfg.t_dist_gamma_shape = str2num(y(2:end));
        case 'BKGD_SCALE'
            cfg.res = str2num(y(2:end));
        case 'USE_SUM_APPROXIMATION'
            cfg.use_sum_approximation = str2num(y(2:end));
        case 'PARAM_T_DIST_TRUNC_LOW'
            cfg.t_dist_trunc_low = str2num(y(2:end));
        case 'PARAM_T_DIST_TRUNC_HIGH'
            cfg.t_dist_trunc_high = str2num(y(2:end));
        case 'BKGD_INTERP_MAX_DIST'
            cfg.B_interp_max_dist = str2num(y(2:end));
        case 'BKGD_CHANGE_THRESH'
            cfg.B_change_th = str2num(y(2:end));
        case 'BKGD_OFFSET'
            cfg.B_offset = str2num(y(2:end));
        case 'BKGD_PARAM_MAX_SUM_THRESH'
            cfg.B_max_sum_th = str2num(y(2:end));
    end
    
    line = fgetl(f);
end

f = fclose(f);


% extract data - B value and segment length for each segment

[b, l] = textread( file, '%d\t%d', 'commentstyle', 'shell');

BL = [max(1,b) l];


end

