% reformat neutral files to matlab binary to be used as memory map objects

neut_dir = '/Users/davidmurphy/GoogleDrive/linked_selection/data/poly/neutral/';
memmap_dir = '/Users/davidmurphy/GoogleDrive/linked_selection/data/memmaps/';



for c=1:22
    
  % downsampled to 1%, thinned for low recombination files:
  % "chr1.thin.neut.rec_0.000_cmmb.random_sample_01"
  input_file = sprintf('%schr%d/chr%d.thin.neut.rec_0.000_cmmb.random_sample_01', neut_dir, c, c);
  output_file = sprintf('%schr%d/chr%d.thin.neut.rec_0.000_cmmb.random_sample_01.dat', memmap_dir, c, c);

  %   input_file = sprintf('%schr%d/chr%d.neutral', neut_dir, c, c); 
%     output_file = sprintf('%schr%d/chr%d.neutral.dat', neut_dir, c, c);

  
  neut_txt = fopen( input_file, 'rt' );

% sample line for neutral data: 16050408	162	15	T	C	1
% read the neutral data file with textscan
% take only the relevant columns 1:3 = 'pos', 'N', 'P'
  neut_data = textscan(neut_txt, '%d%d%d%*[^\n]', 'delimiter','\t');  
  
  fclose(neut_txt);
  
  % write a new matlab binary data file:
  % just write each set of numbers one after the other
  % when creating memory map, we know how many of each type there are,
  % which follows the number of lines in the neutral files
  binary_neut_data = fopen(output_file, 'w');
  
  for r=1:3
    fwrite(binary_neut_data, neut_data{r}, 'int32');
  end 
  
  fclose(binary_neut_data);
  
end

clearvars;
