
function fdata = LS_LoadVariationData( inputfiles, vcfg ) %inputfiles_file

% inputs:
% inputfiles - configuration files containing input file names (one per chromosome)
% outputs:
% fdata      - a struct containing all loaded variation data: genetic maps,
%              polymorphism and local mutation rate estimates


if isstr(inputfiles)
  inputfiles = file2struct( inputfiles );
end
  

%% general inits

% just use this old method for now until sscanf is worked out...
if 0
  base_dir = '/Users/davidmurphy/GoogleDrive/linked_selection/lsm/cluster_mirror/';
else
  base_dir = '/ifs/data/c2b2/gs_lab/dam2214/inf/';
end

chr_features_file = [base_dir 'ch_features/chr_len_all.txt'];
[chr_id, chr_len] = textread(chr_features_file, '%s %d', 'headerlines', 1);

C = length(inputfiles.poly);


for c=1:C
    
    chid = sprintf('chr%d', c);  % assign from cfg_inf
%     chr_features_file = [base_dir sprintf('ch_features/%s_chromosomal_features.txt', chid) ]; % fixed: no longer hard coded, uses chrID from inits
    fdata{c}.chr_id  = chr_id(c);
    fdata{c}.chr_len =  chr_len(c);
    
end




%% load genetic map

for c=1:C
  [genmap{c}.pos, genmap{c}.c, genmap{c}.R]  = textread( inputfiles.genmap{c}, '%d\t%f\t%f', 'commentstyle', 'shell' );
  genmap{c}.file = inputfiles.genmap{c};
  genmap{c}.name = inputfiles.genmap_token;
  
  fdata{c}.genmap = genmap{c};
  fdata{c}.genmaplims = [min(genmap{c}.R) max(genmap{c}.R)] /100;
end



%% load poly data

for c=1:C
    
%   f = fopen( inputfiles.poly{c}, 'rt' );
% 
% % sample line for human data: 16050408	162	15	T	C	1
%   Z = textscan( f, '%d\t%d\t%f\t%s\t%s\t%d', 'commentstyle', '#' );
%     
%   L = length(Z{1});
%   
%   fdata{c}.Poly.pos  = Z{1};
%   fdata{c}.Poly.gpos = applyGenmap2pos( fdata{c}.genmap, fdata{c}.Poly.pos );
  %several index changes here for the correct column numbers: % DM edit
%   idxS    = find( Z{6} == 1 ); % 0=mono, 1=SYN, 2=NS %change idx to 6 since no strand
%   fdata{c}.Poly.idxS    = idxS;
%   fdata{c}.Poly.sample  = Z{2}; %changed
%   fdata{c}.Poly.specS   = [Z{3}(idxS) double(Z{2}(idxS))-Z{3}(idxS)]; %changed
%   fdata{c}.Poly.anc = char(zeros([L 1]));
%   fdata{c}.Poly.der = char(zeros([L 1]));
%   for i=1:L
%   fdata{c}.Poly.anc(i,:)     = Z{4}{i}; %changed
%   fdata{c}.Poly.der(i,:)     = Z{5}{i}; %changed
%   end
% %   fdata{c}.Poly.type    = Z{7};
%   
%   f = fclose(f);

  % REDUCED INPUT: [POS, SAMPLESIZE, SNPCOUNT]
  f = fopen( inputfiles.poly{c}, 'rt' );
  Z = textscan( f, '%d\t%d\t%f', 'commentstyle', '#' );
  L = length(Z{1});
  
  fdata{c}.Poly.pos  = Z{1};
  % gmask = (fdata{c}.Poly.pos > fdata{c}.genmap.pos(1)) & (fdata{c}.Poly.pos <  fdata{c}.genmap.pos(end));
  % fdata{c}.Poly.pos = fdata{c}.Poly.pos(gmask);
  fdata{c}.Poly.gpos = applyGenmap2pos( fdata{c}.genmap, fdata{c}.Poly.pos );
  idxS    = find( Z{3} > 0 );
  fdata{c}.Poly.idxS    = idxS;
  fdata{c}.Poly.sample  = Z{2}; 
  fdata{c}.Poly.specS   = [Z{3}(idxS) double(Z{2}(idxS))-Z{3}(idxS)];   
  
  f = fclose(f);
  
end



%% load mutation rate proxy data

for c=1:C
  f = fopen( inputfiles.mutprox{c}, 'rt' );

% [pos support mutrate]
% OK for my format with no header but remove one column since eyals is 3:
  Z = textscan( f, '%d\t%f', 'commentstyle', '#' ); %remove a float

  
  fdata{c}.MutProx.pos  = Z{1};
  fdata{c}.MutProx.gpos = applyGenmap2pos( fdata{c}.genmap, fdata{c}.MutProx.pos );
  
  % no codons in my mutation rate data:
%   fdata{c}.MutProx.paml_codons  = Z{2}; 
  fdata{c}.MutProx.dS       = Z{2}; %change index to the mutation rate index
 % this looks like it just assigns positions with 0 codon support to be -1
 % just in case it wasnt already done in the datafile:
%   fdata{c}.MutProx.dS( fdata{c}.MutProx.paml_codons<=0 ) = -1;
  
  f = fclose(f);
end

  
end
