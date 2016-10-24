
# LINKED SELECTION INFERENCE: Getting Started
---

0. Basic requirements
1. Introduction
2. Preparing data
3. Configuration File
4. Basic settings
5. Running inference

---
### 1. BASIC REQUIREMENTS
This program is written in MATLAB and requires at least version 20XX or later to run. Additionally, the program calc_bkgd from (McVicker et al., 2009) needs to be compiled on the system running the MATLAB code. This program is written in C and has the following dependencies:

+ gsl-1.16  
+ zlib-1.2.8
+ glib-2.42.2 -- has several dependencies as well:
    - pkg-config 
	- libffi  
	- GNU gettext()  

---
### 2. INTRODUCTION
This program infers the effects of linked selection due to classic sweeps (CS) and background selection (BGS). The models for CS and BGS are detailed in Elyashiv et al., 2015. 

The effects of BGS and CS are estimated genome wide for a range of selection coefficients and across different annotations in a "pre-calculation" step. This step is computationally intensive but for a given selection coefficient and annotation, it is only required once. During the composite likelihood calculation at the core of the program each of these precalculated maps of linked selection are loaded and combined in a weighted average, which is optimized to maximize the composite likelihood. The current configuration outputs the predicted diversity levels on a scale from 0-100 (where 100 = neutral expectation) for both CS, BGS and the combined effects of both.

---
### 3. PREPARING DATA
The minimum required data set to run the inference includes the following:
#### 1. Neutral sites
All neutral sites in the genome based on some criteria. For example, with human data we use sites that have a phastCons 100 vertebrate alignment score of 0, i.e. least conserved (see: [here](http://genome.ucsc.edu/goldenPath/help/phastCons.html)), do NOT fall in any region annotated as coding in the UCSC knownGene genes database, pass 1000 Genomes Project Pilot Mask filters, and have an anotated base in the hg19 reference genome. The reason we check against the 1000 Genomes Project mask file is because we are using 1000 Genomes data for polymorphism data, so we only consider the accessible genome that passes these filters. 

The format for neutral sites data with the current configuration is the following:

    #POS    N   DER
    2731234	162	0
    2731280	162	0
    2731540	162	0
    2731636	162	0
    2731641	162	0
    2731743	162	23
    2731805	162	0
    2731896	162	0
    2732327	162	0
    2732623	162	0

+ POS = position (in 1 based coordinates)
+ N = sample size (i.e. number of chromosomes)
+ DER = derived allele count at the site

#### 2. Genetic maps
A genetic map is required for both BGS and CS as each is a model of the effects of selection due to linkage. Pedigree based maps are preferable to maps inferred from population scale LD because these may be biased by linked selection itself in some cases. The map be formatted such that column 1 indicates position, column 2 indicates the recombination RATE in cM/Mb and column 3 indicates genetic map distance in cM. Note that the relationship between genetic distance and physical position is the following: for a given position, the genetic map distance is the distance from the LAST marker to the current marker.

The format for genetic maps is the following:

    #position	COMBINED_rate(cM/Mb)	Genetic_Map(cM)
    1390135 0.02002002002   0.1048
    1401759 0.0344115622849 0.1052
    1412848 0.928848408333  0.1155
    1423176 0.213013168087  0.1177
    1430916 0.0516795865633 0.1181
    1440776 0.0709939148073 0.1188
    1450323 0.324709332775  0.1219
    1462217 0.378342021187  0.1264
    1471286 0.385930091521  0.1299

#### 3. Annotated segments
Annotated segments are not strictly required in the case of running the inference with CS only, but for BGS a minimum of one segment annotation is required. Examples include "coding", "noncoding", "exonic", "utr", etc. For a given run, segments should not overlap, so if you choose "coding" and "utr" you should not include "exonic" as these segments will overlap. Within a given set of segments, segments must ALSO not overlap -- in general, for a set of indices `i_0, ..., i_(n-1)` for the set of `n` segments, the following conditions must hold:

    segments[i][end] < segments[i+i][start]  (segments do not overlap each other)

and for `i_0, ..., i_n`:

    segments[i][start] < segments[i][end]  (segments have positive length)

Segments are stored in a minimal 3 column BED formatted file such as the following example:

    # chr11 UCSC exon segments from nonredundant coding gene set
    # CHROM    START    END
    chr11	126987	131373
    chr11	131467	131524
    chr11	131752	131920
    chr11	193080	193154
    chr11	193712	193911
    chr11	194418	194573
    chr11	196761	197046
    chr11	197297	197413
    chr11	197561	197761

#### 4. Annotated substitutions
Annotated substitutions are also optional depending on the mode of linked selection you are targetting, but at least one set of annotated substitutions is required for identifying CS. Substitutions are identified using a multi-species alignment or ancestor reconstruction aligned to your target species. For example, in the case of humans the human genome is compared to a human-chimp ancestor based on a 6 species alignment and positions in the human reference genome which differ from the ancestor are called as substitutions. The substitutions are then annotated depending on what kind of annotated segment they fall within. For example, all substitutions which fall within codons and change amino acid translation are annotated as nonsynonymous.

Substitutions file is formatted as in the following example:

        #CHROM 	POS 	DER 	ANC 	ANNO
        chr11	205392	T	C	NS
        chr11	208893	G	A	NS
        chr11	209905	C	A	NS
        chr11	218894	G	C	NS
        chr11	233514	A	G	NS
        chr11	236063	T	C	NS
        chr11	237086	A	G	NS
        chr11	244081	G	C	NS
        chr11	244114	G	A	NS
        chr11	244122	G	A	NS

The fields are:
+ chromosome
+ position
+ derived allele
+ ancestral allele
+ annotation (i.e. NS for "nonsynonymous" in the example)

The only field used with the present configurations is the positions field.

#### 5. Mutation rate estimates
Estimated mutation rates are not required although they improve the inference. Since the we are fitting our model of BGS and CS to observed levels of diversity and diversity levels are proportional to the local mutation rate, we use an estimate of the mutation rate to normalize model predictions of diversity when evaluating the likelihood of the model given the data. Note that, if you use a mutation rate correction, you must also use the correction with the final inference results as well.

We use substitution rates along the genome of an outgroup species (macaque) in windows aligned to our target species (human). For a given window, we estimate the mutation rate by keeping only sites that align with our human neutral sites data set and include an annotated base in both human and macaque. We then count the number of neutral sites where humans have undergone a substitution relative to macaque and divide this by the total number of neutral sites in the window to get the rate. We average the rates over all windows and divide each window rate by this average. For example, if a 10kb window has 4000 neutral sites and 400 substitutions relative to macaque, we would estimate the mutation rate in this window as 0.1. If the average mutation rate for the whole genome is 0.05 then the window used in the example would have a normalized mutation rate of 2.

(For more details on picking the best window size for estimating mutation rates see Elyashiv et al., 2015 supplement 1.6)

Mutation rate file format is as follows:

    # POS	RATE
    180000	1.004114
    200000	0.969781
    220000	1.014950
    240000	1.189234
    260000	1.391958
    280000	1.147105
    300000	1.365735
    320000	1.689727
    360000	1.178376

The position is the _start_ position for the range _(start, start + window)_, so in this case we are using 20KB windows and the rate from 180kb to 200kb is 1.004114.

_NOTE: Even if you do not use a mutation rate correction, you need to specify a constant mutation rate. To use a constant mutation rate, just make a set of files for each chromosome with 2 rows:_
    
    # POS    RATE
    chrom_start    1.0
    chrom_end    1.0

_This will apply a constant rate across your data via interpolation._

#### 6. Chromosome lengths
This simple file contains the name and lengths of each chromosome that will be used for the inference. It should have the following format:

    #chromosome	len
    chr1	249250621
    chr2	243199373
    chr3	198022430
    chr4	191154276
    chr5	180915260
    chr6	171115067
    chr7	159138663
    chr8	146364022
    chr9	141213431

---
### 4. Configuration file
The inference program is designed to run for a range of inputs and settings that can be controlled from a single configuration file. The equivalent of the `main` function calls this configuration script at the beginning of the run to initialize things like file directories, number of annotations, etc. Additionally, there is a default settings file which is meant to be set up _once_ with basic parameters for the species and data set that will be used for the inference, e.g. things such as number of chromosomes, effective population size, range of selection coefficients, etc.

+ The configuration file uses a variable called `base_dir` which is the base directory for all data files. All lower level directories called by the configuration will relate to the `base_dir`. You must assign this variable based on the system you are using to run the inference.
+ One early source of bugs will be misspecified filenames and directories. Go through the configurationFile carefully and make sure that all of the filenames and directories make sense and that they all exist. Some empty directories will be needed for the initial run (they will be filled as you go), so these must be specified as well.
+ Many of the file names are build from smaller units, which are used in different contexts. For example, the genetic map name and map resolution are specified in the following variables:

	genmap_id = 'AA_Map';
	map_res = '10kb';
	
+ Using some of these identifiers, genetic map files are specified:

        genmap_token = sprintf('%s_%s', genmap_id, map_res);
        % write a list of genetic map files
        % specify the full path to the files
        for c=1:C
            map_file = sprintf('maps/%s/%s_%s_window_hg19.gmap',genmap_id, chr_id{c}, genmap_token);
            genmap_files{c} = [base_dir map_file];
        end
	
+ Other files will be constructed using the `genmap_token` variable during the inference.
+ Read through the configuration file and modify it accordingly to suit your data and directory structure. It is heavily commented and most sections should be self-explanatory.
