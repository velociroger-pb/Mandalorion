# Mandalorion #
[![Github release](https://img.shields.io/github/tag/christopher-vollmers/Mandalorion-1.svg?label=Version)](https://github.com/christopher-vollmers/Mandalorion-1/tags)
[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](http://perso.crans.org/besson/LICENSE.html)

## This repository is being heavily modified currently. Stuff should work most of the time but behavior might change somewhat unpredictably. Things should be more stable in September 2021. A stable version of Mandalorion (3.5) lives in Roger Volden's github ##

v3.6.1: This is the Isoform.

Takes R2C2/C3POa or PacBio/ccs/lima data and defines high confidence isoform consensus sequences and alignments. 

## PacBio read conversion/formatting ##

PacBio ccs reads (.fastq) and subreads (.bam) files have to be converted (to .fasta and .fastq, respectively) and formatted to match C3POa naming convention. The bam to fastq conversion can be done with bedtools for example. I do recommend gzipping the resulting fastq files using a multithreaded tool like pigz. Then, you can use the *convertPBcDNAreadsForMandalorion.py* script in utils/ to get ccs reads and subreads into the right naming scheme as well as to subsample subreads to a useable depths. 

## Dependencies ##

- [minimap2](https://github.com/lh3/minimap2)
- [emtrey](https://github.com/rvolden/emtrey) ([go](https://golang.org/dl/))
- [blat source](https://users.soe.ucsc.edu/~kent/src/blatSrc35.zip) or [blat executable](http://hgdownload.soe.ucsc.edu/admin/exe/)
- [medaka](https://github.com/nanoporetech/medaka)

The paths to these will need to be put into your config file [like this.](example_config) If you have the program installed or in your path already, replace the path with the name of the program.

- [mappy](https://pypi.org/project/mappy/)
- [pyabPOA](https://pypi.org/project/pyabpoa/)

After installation with pip3 these should just work

## Usage ##
```bash
python3 Mando.py [OPTIONS]
```

Running with default settings:
```bash
python3 Mando.py -p . -g gencodeV29.gtf -G hg38.fasta -f Consensus_reads.fofn
```

### Options: ###

```bash
                        
  -p PATH, --path PATH  Directory to put output files into
  
  -u UPSTREAM_BUFFER, --upstream_buffer UPSTREAM_BUFFER
                        Defines upstream leniency window for polyA and TSS
                        determination (default 10)
                        
  -d DOWNSTREAM_BUFFER, --downstream_buffer DOWNSTREAM_BUFFER
                        Defines downstream leniency window for polyA and TSS
                        determination (default 50)
                        
  -s SUBSAMPLE_CONSENSUS, --subsample_consensus SUBSAMPLE_CONSENSUS
                        Defines how many random subreads are used to polish
                        isoform consensus sequences (default 500)
                        
  -g GENOME_ANNOTATION, --genome_annotation GENOME_ANNOTATION
                        Genome annotation file (gtf). Is used to identify
                        individual annotated splice sites. If -W is set it is
                        also used to white-list annotated polyA sites
                        
  -G GENOME_SEQUENCE, --genome_sequence GENOME_SEQUENCE
                        Genome file (fasta)
                        
  -r MINIMUM_RATIO, --minimum_ratio MINIMUM_RATIO
                        Proportion of reads that align to a locus required for
                        an isoform (default 0.01)
                        
  -i MINIMUM_INTERNAL_RATIO, --minimum_internal_ratio MINIMUM_INTERNAL_RATIO
                        An isoforms that is completely internal to another isoform 
                        will be discarded unless they reach this ratio compared 
                        to isoform that contains it (default=1)
                        
  -R MINIMUM_READS, --minimum_reads MINIMUM_READS
                        Minimum number of reads to make an isoform (default 5)
                        
  -a ADAPTER_FILE, --adapter_file ADAPTER_FILE
                        Fasta file with 5prime and 3prime adapters. Use if
                        your input reads are not trimmed (default C3POa
                        output) and you want to trimm your isoforms. Don't use
                        if your reads are trimmed, aka default ccs/lima output
                        Will be ignored unless you also use -e to set your
                        sequence ends.
                        
  -f R2C2_CONSENSUS_READS, --R2C2_Consensus_reads R2C2_CONSENSUS_READS
                        Fasta file with R2C2 consensus reads, can be entered
                        as a single file path, a comma separated list of file
                        paths, or a path to a file of filenames file (has to
                        end on .fofn) that contains one file path per line
                        
  -b R2C2_SUBREADS, --R2C2_subreads R2C2_SUBREADS
                        Fastq file(s) with R2C2 subreads, can be entered as a
                        single file path, a comma separated list of file
                        paths, or a path to a file of filenames file (has to
                        end on .fofn) that contains one file path per line. 
                        If Mandalorion is run in consensus mode C, this flag 
                        doesn't have to be used and is ignored if used
                        
  -O OVERHANGS, --overhangs OVERHANGS
                        Defines bounds for unaligned bases on ends. Format:
                        min5prime,max5prime,min3prime,max3prime (default
                        0,40,0,40)
                        
  -t MINIMAP2_THREADS, --minimap2_threads MINIMAP2_THREADS
                        Number of threads to use when running minimap and
                        consensus calling (default 4)
                        
  -e ENDS, --ends ENDS  Ends of your sequences. Use if your input reads are
                        not trimmed and you want to trimm your isoforms, aka
                        default C3POa output. Don't use if your reads are
                        trimmed, aka default ccs/lima output Will be ignored
                        unless you also use -a to set adapter sequences.
                        Format: 5prime,3prime; Smart-seq2 Example: ATGGG,AAAAA
                        
  -I MINIMUM_ISOFORM_LENGTH, --minimum_isoform_length MINIMUM_ISOFORM_LENGTH
                        Minimum length in nt of isoforms that will be
                        considered (default 500)
                        
  -n MINIMUM_FEATURE_COUNT, --minimum_feature_count MINIMUM_FEATURE_COUNT
                        Features (starts,ends, new splice sites) will be
                        considered if they are present in this number of reads
                        (default 2)
                        
  -w SPLICE_SITE_WINDOW, --splice_site_window SPLICE_SITE_WINDOW
                        reads spliced within this number of nucleotides on
                        each side of a splice site will be considered spliced
                        at this site (default 1)
                        
  -A ACUTOFF, --Acutoff ACUTOFF
                        Isoforms with A content of more than this cutoff in a
                        15nt window downstream of their polyA site will be
                        discarded (default 0.5)
                        
  -W, --white_list_polyA
                        If set, polyA sites that fall within +/-20nt of
                        annotated transcript ends will not be filtered
                        regardless of Acutoff set with -A. Annotated
                        transcript ends will be taken from annotation file
                        given with -g
                        
  -S SAM_FILE, --sam_file SAM_FILE
                        If given, Mandalorion will use this file instead of
                        performing its own minimap2 alignment or reads. Careful! If
                        names don't line up between this and the fasta and
                        fastq subread files everything breaks!
                        
  -M MODULES, --Modules MODULES
                        Defines what modules of Mandalorion will be run. By
                        default this includes: A - Alignment, P - .sam to
                        .clean.psl conversion, S - defining splice sites, D -
                        defining isoforms, C - Creating consensus sequences
                        for those isoforms T - Trimming consensus sequences F
                        - Filtering isoforms Q - Quantifying isoforms. Each
                        module need the output of the modules run before it to
                        function properly. Running individual modules can be
                        useful if you want to for example filter with
                        different parameters without rerunning the whole
                        pipeline

  -C CONSENSUSMODE, --consensusMode CONSENSUSMODE
                        Set to P or PC (deault = P). If P, only pyabpoa will
                        be used to make isoform consensus sequences. If PM,
                        medaka will be used as well. PC generates slightly
                        more accurate sequences for R2C2 data but is much slower.
                        Use P is analyzing PacBio data.

  -v, --version         Prints Mandalorion version
```

## Outputs ##

I consider the *Isoforms.filtered.fasta* file the main output of the Mandalorion pipeline. It contains the polished sequences of all isoforms Mandalorion considers very high confidence. Mandalorion also creates *Isoforms.filtered.clean.psl* and *Isoforms.filtered.clean.gtf* files which contain minimap2 alignments of those sequences that had small indels removed. These files are easy to upload to the UCSC Genome Browser to inspect. This version of Mandalorion now also generates a *Isoforms.filtered.clean.quant* file which contains the number of R2C2/PacBio reads that associate which each isoform for each given fasta files.

## Utils ##
These are the scripts used to do haplotype phasing and HLA analysis.
