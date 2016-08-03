

#include <stdio.h>
#include <string.h>
#include <glib.h>

#include "bkgd_files.h"


long	count_relevant_lines(FILE *f, int verify_name, char *chr_name)
{
	char	line[1024], name[1024];
	long curpos, valid_lines = 0, relevant_lines = 0, t1, t2;
	
	curpos = ftell(f);
	
	fseek(f, 0, SEEK_SET);

	while ( fgets ( line, sizeof(line), f ) != NULL ) 
		if(line[0]!='\n' && line[0]!='#')
        {
			valid_lines++;
            if (verify_name)
            {
               sscanf(line, "%s\t%ld\t%ld", name, &t1, &t2);
               if( strcmp(name, chr_name) == 0 )
                   relevant_lines++;
            }
            else
                relevant_lines++;
        }
    
	fseek(f, curpos, SEEK_SET);

    if (verify_name)
       fprintf(stderr, "reading list of conserved region: read %ld regions associated with chromsome %s; ignored %ld regions associated with other chromosomes.\n", relevant_lines, chr_name, (valid_lines-relevant_lines));
	return relevant_lines;
}

SeqFeature* load_conserved_from_file(char *filename, char *chr_name, long *pn_sf)
{
	FILE *f_in;
	SeqFeature* sf;
	long	i=0, t1, t2;
	char	line[1024], name[1024];

	f_in = fopen(filename, "rt");

	/* counting the number of valid lines = features */
	*pn_sf = count_relevant_lines(f_in, 1, chr_name);

	/* allocate memory for features */
	sf = g_new(SeqFeature, *pn_sf);

	while ( fgets ( line, sizeof(line), f_in ) != NULL ) 
		if(line[0]!= '\n' && line[0]!= '#'){
			seqfeat_init(&sf[i]);
			sscanf(line, "%s\t%ld\t%ld", name, &t1, &t2); // Guy: changed d to ld
			if( strcmp(name, chr_name) == 0 )
            {
               sf[i].c.start = t1;
               sf[i].c.end   = t2;
               i++;
            }
		}

	fclose(f_in);

	return sf;
}



SeqFeature* load_genetic_map_from_file(char *filename, Chromosome *chr, long *pn_sf)
{
	FILE *f_in;
	SeqFeature* sf;
	long	start, i=0;
	char	line[1024];
	float	rate, cmorgans;
	long bigo;

	f_in = fopen(filename, "rt");

	/* counting the number of valid lines = features */
	*pn_sf = count_relevant_lines(f_in, 0, NULL);

	/* allocate memory for features */
	sf = g_new(SeqFeature, *pn_sf);

//	// read headers line
//	fgets ( line, sizeof(line), f_in );

	while ( fgets ( line, sizeof(line), f_in ) != NULL ) 
		if(line[0]!= '\n' && line[0]!= '#'){
			seqfeat_init(&sf[i]);
			bigo = sscanf(line, "%ld %f %f", &start, &rate, &cmorgans); // Guy : changed d to ld
			sf[i].c.start = start;
			if(i>1)	sf[i-1].c.end = start-1;
			sf[i].score = rate;
			i+=1;
		}
	
	sf[i-1].c.end = chr->len;

	fclose(f_in);

	return sf;
}



/**
 * Reads a specific chromosome's features from a file
 * GUY: removed this one
long	get_chr_features(Chromosome* chr, char* chr_token, char* chr_features_file)
{
	FILE *f_in;
	char	line[1024], sztemp[100];

	f_in = fopen(chr_features_file, "rt");

	while ( fgets ( line, sizeof(line), f_in ) != NULL ) 
		if(line[0]!= '\n' && line[0]!= '#'){
			sscanf(line, "%s", sztemp);
			if(strcmp(sztemp, chr_token)==0)
				break;
		}

	fclose(f_in);

	chr->name = g_new(char, 256);
	sscanf(line, "%s\t%ld", chr->name, &chr->len); // Guy : changed d to ld

	chr->id = -1;
	chr->assembly = NULL;

	return 0;
}
 */

long	print_params(FILE *f_out, Chromosome *chr, char *params_file)
{
	FILE *f_params;
	char	line[1024];


	fprintf(f_out, "#%s\t%ld\n", chr->name, chr->len); // Guy : changed d to ld

	f_params = fopen(params_file, "rt");

	while ( fgets ( line, sizeof(line), f_params ) != NULL ) 
		if(line[0]!= '\n' && line[0]!= '#')
			fprintf(f_out, "#%s", line);

	fclose(f_params);

	return 0;
}
