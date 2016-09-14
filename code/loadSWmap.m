
function [SW, SWfake, pGrid, cfg] = loadSWmap( file, cfg )

% extract parameters from the commented lines

f = fopen(file,'rt');

if f==-1
  error('SW-map file not found');
end

line = fgetl(f);

[cfg.chr_id, cfg.chr_len] = sscanf(line, '%s %d');

  
while (~isempty(line) & (line~=-1) & (line(1)=='#'))
    [x, y] = strtok(line(2:end), '=');
    switch x
        case 'OUTPUT_FILE'
          ofile = y(2:end);
          i1 = find(ofile=='\',1,'last');
          if ~isempty(i1)
            cfg.output_file = ofile(i1+1:end);
            cfg.output_dir = ofile(1:i1-1);
          else
            cfg.output_file = ofile;
            cfg.output_dir = '';
          end
        case 'OUTPUT_TOKEN'
            cfg.output_token = y(2:end);
        case 'CHROMOSOME_NAME'
            cfg.chr_id = y(2:end);
        case 'CHROMOSOME_LENGTH'
            cfg.chr_len = str2num(y(2:end));
        case 'RECOMB_RATE_TABLE'
            cfg.rec_table = y(2:end);
        case 'RECOMB_RATE_SCALE'
            cfg.rec_table_scale = str2num(y(2:end));
        case 'FOCALS_TABLE'
            cfg.focals_table = y(2:end);
        case 'PARAM_NE0'
            cfg.Ne0 = str2num(y(2:end));
        case 'PARAM_S'
            cfg.S = str2num(y(2:end));
        case 'PARAM_S_DIST_TYPE'
            cfg.s_dist_type = y(2:end);
        case 'TTF_APPROXIMATION'
            cfg.trap_aprx = y(2:end);
        case 'USE_FAKE_SUBSTITUTIONS'
            cfg.use_fake_subs = str2num(y(2:end));
        case 'StopSum'
            cfg.StopSum = str2num(y(2:end));
        case 'INTERP_METHOD'
            cfg.InterpMethod = y(2:end);
        case 'MAX_DIST_SCALED'
            cfg.gMaxDistScaled = str2num(y(2:end));
        case 'MAX_DIST'
            cfg.gMaxDist = str2num(y(2:end));
    end
    
    line = fgetl(f);
end

f = fclose(f);


f = fopen( file, 'rt' );
% Z = textscan( f, '%d\t%f\t%f', 'commentstyle', {'#'}); %'shell'
% no 'fake' field:
Z = textscan( f, '%d\t%f', 'commentstyle', '#');
f = fclose(f);

pGrid = Z{1};  SW = double(Z{2});  SWfake = zeros([length(pGrid) 1]);

% [pGrid, SW, SWfake] = textread( file, '%d\t%.3f\t%.3f', 'commentstyle', 'shell');


end

