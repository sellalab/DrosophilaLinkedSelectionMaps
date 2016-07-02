
function s = file2struct( filename )

% currently supported fields: numeric matrices, (single) strings, string matrices, numeric-matrix matrices.
% unsupported: struct arrays, mixed matrices/strings cell arrays.
% syntax: ';' after the field name separate it from its value, if it's a single entry rather than a 1X1 cell matrix.
%         '|' after the field name separate it from its value for all cell matrices.
%         ',' separates matrix entries.
%         ';' separates matrix rows.
%         '!' separates cell matrix entries.
%         '|' separates cell matrix rows.
%         'UNSUPPORTED' = unsupported field type, a place-holder
%
% NOTE! you cannot have a 1X1 cell array - it will be considered as a single number/string/numeric-matrix
%
% example:
% .cfg.rec_th_L; 7.500000e-01                                % floating point number
% .cfg.rec_th_H; Inf                                         % Infinite (there's also 'NaN')
% .cfg.inf.init_params;                                      % and empty array
% .cfg.inf.fixed_params; 1, 1, 0; 1, 1, 1; 1, 1, 1; 1, 1, 0  % integers 2D matrix
% .cfg.GEs.CalcSW.InterpMethod; linear                       % a string
% .cfg.chr_ids; 2L | 2R | 3L | 3R | X                        % 4X1 matrix of strings
% .params.SWparams; UNSUPPORTED                              % unsupported field type
% .params.BS.w_t; 0, 1: 0, 1: 0, 1: 0, 1                     % 1X2 cell matrix of matrices


s = [];

f = fopen( filename, 'rt' );


lien = fgetl( f );

while isempty(lien) | (~isempty(lien) & lien~= -1)
  
  if ~isempty(lien) & ~strcmp(lien,'')
    
  lien = lien(lien~=' ');
  
  i = find( lien==';' | lien=='|' , 1);
  assert( ~isnan(i), 'invalid syntax - no proper field-value separator' );
  
  if lien(i)==';'
    is1x1cell = 0;
  else
    is1x1cell = 1;
  end
  
  %   if isempty(i)
  %
  %     values = [];
  %     prefix = lien;
  %     subfields = strsplit( lien, {'.'} );
  %
  %   else
  
  prefix = lien(1:i-1);
  subfields = strsplito( lien(1:i-1), '.' );
  %   subfields = strsplit( lien(1:i-1), {'.'}, 'CollapseDelimiters', false   );
  
  % if strcmp(prefix,'.inf.predef_idx_train')
  %   bigo = 7;
  % end

  cell_rows = strsplito( lien(i+1:end), '|' );
  %   cell_rows = strsplit( lien(i+1:end), {'|'}, 'CollapseDelimiters', false  );
  
  vvalues = [];
  for l=1:length(cell_rows)
    
    cells = strsplito( cell_rows{l}, '!' );
    %     cells = strsplit( cell_rows{l}, {'!'}, 'CollapseDelimiters', false );
    
    for k=1:length(cells)
      
      values = strsplito( cells{k}, ',;' );
      %       values = strsplit( cells{k}, {',',';'}, 'CollapseDelimiters', false   );
      h = 1 + sum( cells{k}==';' );
      w = length( values )/h;
      
      nums = cellfun( @str2num, values, 'UniformOutput', false );
      
      
      
      if sum( cellfun( @isempty,nums ) ) == 0
        % nums
        values = cell2mat( nums );
      elseif length(values)==1
        if isstr(values{1}) & strcmp(values{1}, 'UNSUPPORTED')
          values = [];
        else
          values = values{1};
        end
      end
      
      if h~=1 | w~=1
        values = reshape( values, [w h] )'; % this transpose-patch is needed because the series of numbers is read as a 1Xn vector
      end
      %   end
      
      if strcmp(cells{k},'')
        values = [];
      end
      vvalues{l,k} = values;
      
    end
  end
  % if the cell matrix has a single cell and its contents is a string, remove cell encapsulation
  if prod(size(vvalues))==1 & is1x1cell==0
    vvalues = vvalues{1,1};
  end
  
  s = setfield_rec( s, subfields{2}, subfields(3:end), vvalues );
  
  end
  lien = fgetl( f );
  
end

f = fclose( f );

end

%%

function s = setfield_rec(s, field, subfields, value)

% s = [];

if length(subfields)==0
  s = setfield(s, field, value);
else
  if isfield(s, field)
    sub_s = getfield(s,field);
  else
    sub_s = [];
  end
  s = setfield(s, field, setfield_rec(sub_s, subfields{1}, subfields(2:end), value) );
end

end



%%

function splitted = strsplito( str, delimiters )

if isempty(str)
  splitted = {''};
else
  splitted = textscan( str, '%s', 'delimiter', delimiters, 'MultipleDelimsAsOne', 0 );
  splitted = splitted{1};
end

end

