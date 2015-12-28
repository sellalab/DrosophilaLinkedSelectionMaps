
function Annots   = LS_LoadSelectionAnnotations(...
  chr_features_file,...
  chr_id,...
  genmap_files, genmap_token,...
  SW_anno_files, SW_anno_tokens,...
  BS_anno_files, BS_anno_tokens )




%% general inits

[ff_chr_id, ff_chr_len] = textread(chr_features_file, '%s\t%d', 'commentstyle', 'shell' ); %, 'headerlines', 1

C = length(chr_id);

for c=1:C
  ichr(c) = find(strcmp( ff_chr_id, chr_id{c} ));
  chr_len(c) = ff_chr_len(ichr(c));
end


% % load genetic map
%
% for c=1:C
%   [genmap{c}.pos, genmap{c}.c, genmap{c}.R]  = textread( genmap_files{c}, '%d\t%f\t%f', 'headerlines', 1 );
%   genmap{c}.file = genmap_files{c};
%   genmap{c}.name = genmap_token;
% end



%% SW annotations

if ~isempty(SW_anno_files)
  for c=1:C
    for b=1:length(SW_anno_files{c})
      if ~isempty(SW_anno_files{c}{b})
        Annots.SW{c}{b}.output_token  = SW_anno_tokens{b};
        Annots.SW{c}{b}.name         = [chr_id{c} '_' SW_anno_tokens{b}];
        Annots.SW{c}{b}.file  = SW_anno_files{c}{b};
        
%         [Annots.SW{c}{b}.focals.pos, Annots.SW{c}{b}.focals.strand, Annots.SW{c}{b}.focals.isfake] = ...
%           textread( Annots.SW{c}{b}.file, '%d\t%c\t%d', 'commentstyle', 'shell' );


        % edited 12/27: we don't care about strand or "isfake"
        % formatting of my files: chr1	874672	G	C	NS
        [~, Annots.SW{c}{b}.focals.pos, ~, ~, ~] = ...
          textread( Annots.SW{c}{b}.file, '%s\t%d\t%s\t%s\t%s', 'commentstyle', 'shell' );
        
        %       Annots.SW{c}{b}.focals.gpos  = applyGenmap2pos( genmap{c}, Annots.SW{c}{b}.focals.pos );
        Annots.SW{c}{b}.anno_len           = length( Annots.SW{c}{b}.focals.pos ); 
        % edit 12/27: just fake the other fields above for now, they dont 
        % really matter:
        for i=1:Annots.SW{c}{b}.anno_len
            Annots.SW{c}{b}.focals.strand{i} = '+';
            Annots.SW{c}{b}.focals.isfake(i) = 0;  % none of the subs are fake
        end
        
      end
    end
  end
end




%% BS annotations
if ~isempty(BS_anno_files)
  for c=1:C
    for b=1:length(BS_anno_tokens)
      if ~isempty(BS_anno_files{c}{b})
        Annots.BS{c}{b}.output_token  = BS_anno_tokens{b};
        Annots.BS{c}{b}.name         = [chr_id{c} '_' BS_anno_tokens{b}];
        Annots.BS{c}{b}.file  = BS_anno_files{c}{b};
        
        [~, Annots.BS{c}{b}.istart, Annots.BS{c}{b}.iend] = textread( Annots.BS{c}{b}.file, '%s\t%d\t%d' );
        
        %       Annots.BS{c}{b}.gistart  = applyGenmap2pos( genmap{c}, Annots.BS{c}{b}.istart );
        %       Annots.BS{c}{b}.giend    = applyGenmap2pos( genmap{c}, Annots.BS{c}{b}.iend );
        Annots.BS{c}{b}.anno_len           = sum( Annots.BS{c}{b}.iend-Annots.BS{c}{b}.istart+1 );
      end
    end
  end
end




end

