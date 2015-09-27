
function struct2file( s, filename )

% currently supported fields: numeric matrices, (single) strings, string matrices, numeric-matrix matrices.
% unsupported: struct arrays, mixed matrices/strings cell arrays.
% syntax: ';' after the field name separate it from its value, if it's a single entry rather than a 1X1 cell matrix.
%         '|' after the field name separate it from its value for all cell matrices.
%         ',' separates matrix entries.
%         ';' separates matrix rows.
%         '!' separates cell matrix entries.
%         '|' separates cell matrix rows.
%         'UNSUPPORTED' = unsupported field type, a place-holder

f = fopen( filename, 'wt' );

struct2file_rec( s, f, '' );
% struct2file_rec( s, f, inputname(1) );

f = fclose( f );

end


function struct2file_rec( s, f, prefix )

fields = fieldnames(s)';

mstruct = structfun( @isstruct, s );


for i=1:length(mstruct)
  
  cs = getfield(s, fields{i});
  
  if mstruct(i)
    struct2file_rec( cs, f, [prefix '.' fields{i}] );
  else
    % allow the following fields: numerical 1D/2D array, character array, string, string cell array
    
    dims = length(size( cs ));
    assert( dims <= 2, 'unsupported field dimensionality');
%     [h w] = size( cs );
    
    % print field full name
    fprintf( f, [prefix '.' fields{i}] );
%     fprintf( f, [prefix '.' fields{i} ';'] );
    
    if isnumeric( cs ) % numeric matrix
      fprintf( f, ';' );
      for n=1:size(cs,1)
        if n>1, fprintf( f, ' ;');, end
        if (cs(n,1)==round(cs(n,1))),   fprintf( f, ' %d', cs(n,1) );,   else, fprintf( f, ' %e', cs(n,1) );,   end
        for m=2:size(cs,2)
          if (cs(n,m)==round(cs(n,m))), fprintf( f, ' , %d', cs(n,m) );, else, fprintf( f, ' , %e', cs(n,m) );, end
        end
      end
    elseif iscellstr( cs ) % string cell matrix
      fprintf( f, '|' );
      for n=1:size(cs,1)
        if n>1, fprintf( f, ' |');, end
        fprintf( f, ' %s', cs{n,1} );
        for m=2:size(cs,2)
          fprintf( f, ' ! %s', cs{n,m} );
        end
      end
    elseif ischar( cs ) % single string
      fprintf( f, ';' );
      fprintf( f, ' %s', cs );
    else
      ismatrixcellarray = iscell(cs);
      for k=1:length(cs)
        ismatrixcellarray = ismatrixcellarray & isnumeric(cs{k});
      end
      if ismatrixcellarray % cell matrix of numeric matrices
        fprintf( f, '|' );
        for k=1:size(cs,1)
          if k>1, fprintf( f, ' |');, end
          for l=1:size(cs,2)
            %           [h w] = size( cs{k} );
            for n=1:size(cs{k},1) % h
              if n>1, fprintf( f, ' ;');, end
              if (cs{k}(n,1)==round(cs{k}(n,1))),   fprintf( f, ' %d', cs{k}(n,1) );,   else, fprintf( f, ' %e', cs{k}(n,1) );,   end
              for m=2:size(cs{k},2) % w
                if (cs{k}(n,m)==round(cs{k}(n,m))), fprintf( f, ' , %d', cs{k}(n,m) );, else, fprintf( f, ' , %e', cs{k}(n,m) );, end
              end
            end
            if l<size(cs,2), fprintf( f, ' !' );, end
          end
        end
      else
        fprintf( f, '; UNSUPPORTED' );
      end
    end
    
    fprintf( f, '\n' );
  end
  
end


end


%%
% function struct2File( s, fileName, varargin )
% %Write struct to text file. The data in the struct can be both numbers and
% %strings. The first line in the file will be a row with the headers/names
% %of the columns.
% %Ex: struct2File( s, 'c:/test.txt' );
% %Ex: struct2File( s, 'c:/test.txt', 'precision', 3 ); %specify precision
% %Ex: struct2File( s, 'c:/test.txt', 'promptOverWrite', false );
%
%
% [varargin,align]=getarg(varargin,'align',false);
% align=align{:};
% [varargin,delimiter]=getarg(varargin,'delimiter','\t');
% delimiter=delimiter{:};
% [varargin,units]=getarg(varargin,'units','');
% units=units{:};
% [varargin,promptOverWrite]=getarg(varargin,'promptOverWrite',false);
% promptOverWrite=promptOverWrite{:};
% [varargin,precision]=getarg(varargin,'precision',6);
% precision=precision{:};
% [varargin,Sort]=getarg(varargin,'sort',true);
% Sort=Sort{:};
%
% if ~isempty(varargin)
%     error('Unknown optional arguments specified');
% end
%
% fields = fieldnames(s)';
%
% if ~isempty(units)
%     if numel(units)~=numel(fields)
%         error('The number of units specified doesn not match the number of fields in the struct');
%     end
% end
%
% if exist(fileName,'file')==2 && promptOverWrite
%     res = questdlg('File exists, overwrite?','', ...
%         'Yes', 'No', 'Yes');
%     if strcmpi(res,'No')
%         disp('Aborted');
%         return;
%     end
% end
%
% data=cell(numel(s),numel(fields));
% for k=1:numel(fields)
%     fn=fields{k};
%     data(:,k) = {s.(fn)};
% end
% if size(units,2)==1
%     units = units';
% end
% if ~isempty(units)
%     data=[units;data];
% end
% data=[fields;data];
%
% if Sort
%     [fields,ind] = sort(fields);
%     data = data(:,ind);
% end
%
% %ex1 = {'a' 1 12 123; 'ab' 4 5 6; 'abc' 7 8 9};
% ex_func3 = @(input)ex_func(input,precision);
% ex2 = cellfun(ex_func3,data,'UniformOutput',0);
% if align
%     size_ex2 = cellfun(@length,ex2,'UniformOutput',0);
%     str_length = max(max(cell2mat(size_ex2)))+1;
%     ex2 = cellfun(@(x) ex_func2(x,str_length),ex2,'uniformoutput',0);
%     ex2 = cell2mat(ex2);
% end
%
% fid = fopen(fileName,'wt');
% if fid==-1
%     error('Could not open %s');
% end
%
% if iscell(ex2)
%     [m,n]=size(ex2);
%     for i=1:m
%         for j=1:n
%             fprintf(fid,'%s',ex2{i,j});
%             if j<n
%                 fprintf(fid,delimiter);
%             end
%         end
%         if i<m
%             fprintf(fid,'\n');
%         end
%     end
% else
%     m=size(ex2,1);
%     for i=1:m
%         fprintf(fid,'%s',strtrim(ex2(i,:)));
%         if i<m
%             fprintf(fid,'\n');
%         end
%     end
% end
% fclose(fid);
%
%
%
%
% function [ out ] = ex_func( in, prec )
%
% if iscell(in)
%     in = in{:};
% end
%
% in_datatype = class(in);
%
% switch in_datatype
%     case 'char'
%         out = in;
%     case 'double'
%         out = num2str(in, prec);
%     otherwise
%         error('Unknown type');
% end
%
% function [ out ] = ex_func2( in, str_length )
%
% a = length(in);
% out = [char(32*ones(1,str_length-a)), in];
