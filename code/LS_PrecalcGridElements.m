function [GEs, Annots]   = LS_PrecalcGridElements( inoutfiles, cfg )

% inputs:
% inoutfiles - a configuration file holding input and output file names
% cfg        - a configuration file holding parameters for the creation of grid elements
% outputs:
% GEs - a struct containing calculated grid elements
% Annots - a struct contaning annotation lists (just for debug - not needed for further processing)



if isstr(inoutfiles)  % the build GE files
  inoutfiles = file2struct(inoutfiles);
end
if isstr(cfg)  % the global config parameters, held in the config for the whole run
  cfg = file2struct(cfg);
end



outdir            = inoutfiles.outdir;  % where the grid files are placed
outfile_pref      = inoutfiles.outfile_pref;
chr_features_file = inoutfiles.chr_features_file;
chr_id            = inoutfiles.chr_id;
pos_grid_files    = inoutfiles.pos_grid_files;
genmap_files      = inoutfiles.genmap_files;
genmap_token      = inoutfiles.genmap_token;

if isfield(inoutfiles,'SW_anno_tokens'), SW_anno_tokens    = inoutfiles.SW_anno_tokens;, else, SW_anno_tokens = [];, end
if isfield(inoutfiles,'BS_anno_tokens'), BS_anno_tokens    = inoutfiles.BS_anno_tokens;, else, BS_anno_tokens = [];, end

SW_anno_files = [];
BS_anno_files = [];
C = length(chr_id);
for c=1:C % DM edit: only 2 annos to run currently:
  if isfield(inoutfiles, 'SW_anno1_files'),  SW_anno_files{c}{1} = inoutfiles.SW_anno1_files{c};,  end 
  if isfield(inoutfiles, 'SW_anno2_files'),  SW_anno_files{c}{2} = inoutfiles.SW_anno2_files{c};,  end
  if isfield(inoutfiles, 'SW_anno3_files'),  SW_anno_files{c}{3} = inoutfiles.SW_anno3_files{c};,  end
  if isfield(inoutfiles, 'SW_anno4_files'),  SW_anno_files{c}{4} = inoutfiles.SW_anno4_files{c};,  end
end

for c=1:C
  if isfield(inoutfiles, 'BS_anno1_files'),  BS_anno_files{c}{1} = inoutfiles.BS_anno1_files{c};,  end 
  if isfield(inoutfiles, 'BS_anno2_files'),  BS_anno_files{c}{2} = inoutfiles.BS_anno2_files{c};,  end
  if isfield(inoutfiles, 'BS_anno3_files'),  BS_anno_files{c}{3} = inoutfiles.BS_anno3_files{c};,  end
  if isfield(inoutfiles, 'BS_anno4_files'),  BS_anno_files{c}{4} = inoutfiles.BS_anno4_files{c};,  end
  if isfield(inoutfiles, 'BS_anno4_files'),  BS_anno_files{c}{5} = inoutfiles.BS_anno5_files{c};,  end

end




%% general inits

GEs.cfg = cfg.GEs;

[ff_chr_id, ff_chr_len] = textread( chr_features_file, '%s %d', 'commentstyle', 'shell' ); %, 'headerlines', 1

% C = length(chr_id);

for c=1:C
  ichr(c) = find(strcmp( ff_chr_id, chr_id{c} ));
  chr_len(c) = ff_chr_len(ichr(c));
end


% load genetic map

for c=1:C
  [genmap{c}.pos, genmap{c}.c, genmap{c}.R]  = textread( genmap_files{c}, '%d\t%f\t%f', 'commentstyle', 'shell'); % , 'headerlines', 1 DM edit for my format
  genmap{c}.file = genmap_files{c};
  genmap{c}.name = genmap_token;
end



%% load selection annotations (SW & BS)

% NOTE: the annotations are loaded to 'Annots' in parallel to their
% usage as inputs to 'generateBbase' through the files
Annots   = LS_LoadSelectionAnnotations(...  % a large struct containing all of the ns and nc subs and many chromosomal features
  chr_features_file,...
  chr_id,...
  genmap_files, genmap_token,...
  SW_anno_files, SW_anno_tokens,...
  BS_anno_files, BS_anno_tokens );


%% SW grid elements

if ~isempty(SW_anno_files)
  CalcSW = GEs.cfg.CalcSW;
  
  % will eval true with default settings:
  if ~CalcSW.skip_generate_maps
%     'creating complete SW maps'
    % calculate and save SW grid elements
    
    for c=1:C
      g1(c) = min(genmap{c}.R)/100; % map limits in Morgans
      g2(c) = max(genmap{c}.R)/100;
    end
    % DM note: the min_mapdist is MUCH smaller than any of the s / 100 values

    gdeltaNF = max( CalcSW.FE_grid/100, CalcSW.min_mapdist_SW_spat_grid ); %3*(sum(g2)-sum(g1))/sum(chr_len)
    
    for c=1:C
      
      for k=1:length(gdeltaNF)
        
        % read positions from the neutral sites files:
        f = fopen( pos_grid_files{c}, 'rt' );
        
        % edit 12/27: using neutral files as pos_grid files:
        % formatting: 780428	162	0
        Z = textscan( f, '%d\t%d\t%d', 'commentstyle', 'shell' );
        pos  = Z{1};
          
%         pos = textread( pos_grid_files{c}, '%d', 'commentstyle', 'shell' ); % reads the SW_bases file -- change to read from main poly file...
        gpos = applyGenmap2pos( genmap{c}, pos ); % convert all distances to genetic units (in Morgans)
        
        [gFocGrid{c}.gpos{k}, idx] = thinnedToGrid( gpos, [gdeltaNF(k) g1(c)-gdeltaNF(k) g2(c)+gdeltaNF(k)] );
        %   gFocGrid{c}.gpos{k} = [g1(c)-gdeltaNF(k):gdeltaNF(k):g2(c)+gdeltaNF(k)];
        
        gFocGrid{c}.pos{k}  = pos(idx);
        %     gFocGrid{c}.pos{k}  = max(1, min(chr_len(c), round(applyGenmap2pos( genmap{c}, gFocGrid{c}.gpos{k}, 1 )) ));
        
      end
    end
  end
 
  
  for c=1:C
    
    CalcSW.chr_id      = chr_id{c};
    CalcSW.chr_len     = chr_len(c);
    %     CalcSW.pos_grid_file = pos_grid_files{c};
    CalcSW.rec_table = genmap_files{c};
    CalcSW.rec_table_scale = 10^-8;
    
    for b=1:length(SW_anno_files{c})
      
      CalcSW.output_token = Annots.SW{c}{b}.output_token; %SW_anno_tokens{b};
      CalcSW.focals_table = Annots.SW{c}{b}.file; %SW_anno_files{c}{b};
      CalcSW.name         = [genmap{c}.name '_' Annots.SW{c}{b}.name];
      grid_file_pref = sprintf('%s/%s%s', outdir, outfile_pref, CalcSW.name);
      
%       Annots.SW{c}{b}.file                                        = CalcSW.focals_table;
%       [Annots.SW{c}{b}.focals.pos, Annots.SW{c}{b}.focals.strand, Annots.SW{c}{b}.focals.isfake] = ...
%         textread( Annots.SW{c}{b}.file, '%d\t%c\t%d', 'commentstyle', 'shell' );
%       Annots.SW{c}{b}.focals.gpos  = applyGenmap2pos( genmap{c}, Annots.SW{c}{b}.focals.pos );
      
      if ~CalcSW.skip_generate_maps % make the maps:
        
        Annots.SW{c}{b}.focals.gpos = applyGenmap2pos( genmap{c}, Annots.SW{c}{b}.focals.pos );
        % the sweep coefficients will constitute the lookup table of values
        % used in the composite likelihood
        
        SWbase{c,b}      = SwCoef( CalcSW.FE_grid, CalcSW.Ne0, gFocGrid{c}, {Annots.SW{c}{b}.focals.gpos(Annots.SW{c}{b}.focals.isfake==0)}, CalcSW );
        SWbase_fake      = SwCoef( CalcSW.FE_grid, CalcSW.Ne0, gFocGrid{c}, {Annots.SW{c}{b}.focals.gpos(Annots.SW{c}{b}.focals.isfake==1)}, CalcSW ); % DM edit-- switch off since there are no fake focals right now...

        SWbase{c,b}.gSWj_fake = SWbase_fake.gSWj;
               
        te = SaveSWBase( grid_file_pref, gFocGrid{c}, SWbase{c,b} );
        CalcSW.grid_files = te.outputfiles;
        
      else % load the maps from files:
        
        [SWbase{c,b}, gFocGrid{c}]      = LoadSWBase( grid_file_pref, CalcSW );
        for k=1:length(gFocGrid{c}.pos)
          gFocGrid{c}.gpos{k} = applyGenmap2pos( genmap{c}, gFocGrid{c}.pos{k} );
        end
        
      end
      
    end
  end

    
  GEs.gFocGrid   = gFocGrid;
  GEs.SWbase     = SWbase;
  if CalcSW.skip_generate_maps
    GEs.cfg.CalcSW  = CalcSW; % we assume here that the cfg struct is the same for all elements (in terms of parameters, not underlying files and annotations)
  end
  
end





%% BS grid elements
if ~isempty(BS_anno_files)
    CalcBS = GEs.cfg.CalcBS;

    % prepare BGS configuration
    CalcBS.chr_features       = chr_features_file;
    CalcBS.rec_table_scale    = 10^-8;
    CalcBS.output_dir         = [outdir '/'];
    CalcBS.output_pref        = outfile_pref;

    for c=1:C
        CalcBS.rec_table          = genmap_files{c};
        CalcBS.chr_id      = chr_id{c};
        CalcBS.chr_len     = chr_len(c);

        for b=1:length(BS_anno_tokens) 

          CalcBS.output_token      = Annots.BS{c}{b}.output_token; %BS_anno_tokens{b};
          % is this right? cons_table is pointing to coding and noncoding segments...
          CalcBS.cons_table   = Annots.BS{c}{b}.file; %BS_anno_files{c}{b}
          CalcBS.name         = [genmap{c}.name '_' Annots.BS{c}{b}.name];

        %       % NOTE: the annotations are loaded to 'Annots' in parallel to their
        %       % usage as inputs to 'generateBbase' through the files
        %       [~, Annots.BS{c}{b}.istart, Annots.BS{c}{b}.iend] = textread(CalcBS.cons_table, '%s\t%d\t%d'); %Annots.BS{c}{b}.file
        %       Annots.BS{c}{b}.anno_len           = sum(Annots.BS{c}{b}.iend-Annots.BS{c}{b}.istart+1);

          % this function EITHER LOADS the maps from files OR CALCULATES&WRITES them to files
          % important edit 12/17 -- now CalcBS.FE_grid is loaded CalcBS.FE_grid{b} because each anno has its own grid
          % [CalcBS, BSbase{c,b}]  = generateBbase(CalcBS, CalcBS.FE_grid{b}, CalcBS.skip_generate_maps);
          [CalcBS, BSbase{c,b}]  = generateBbase(CalcBS, CalcBS.FE_grid(b,:), CalcBS.skip_generate_maps);
          BSbase{c,b}.cfg = CalcBS;
          BSbase{c,b}.cfg.anno_len = Annots.BS{c}{b}.anno_len;

        end
    end


        GEs.BSbase      = BSbase;
    if CalcBS.skip_generate_maps
        GEs.cfg.CalcBS  = CalcBS; % we assume here that the cfg struct is the same for all elements (in terms of parameters, not underlying files and annotations)
    end
end




end

