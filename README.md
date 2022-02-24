# Mandalorion #
[![Github release](https://img.shields.io/github/tag/christopher-vollmers/Mandalorion-1.svg?label=Version)](https://github.com/christopher-vollmers/Mandalorion-1/tags)
[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](http://perso.crans.org/besson/LICENSE.html)


v3.6.2: This is the Isoform.

Takes R2C2/C3POa or PacBio/ccs/lima data and defines high confidence isoform consensus sequences and alignments. 

## PacBio read conversion/formatting ##

PacBio ccs reads (.fastq) have to be converted to fasta files and formatted to match C3POa naming convention. 
Use the *convertPBcDNAreadsForMandalorion.py* script in utils/ for this purpose.

## PolyA tail removal

It is best to run Mandalorion with reads that have no adapter sequences (trimmed), no polyA tails, and in 5'->3' orientation. By default, C3POa (v2.3.0) does trimm reads but doesn't do polyA removal and, depending on the adapters used, may not put reads in 5'->3' orientation..   
Therefore removePolyA.py script in utils/ should be used to get both PacBio (depending on preprocessing) and R2C2 reads ready. 
The script first removes a fixed number of bases (can't be zero) from the ends of the read and then removes the polyA tail on the 3' end of the read. Depending on the polyA position it will also reorient the read in 5'->3' direction. 

## Dependencies ##

- [minimap2](https://github.com/lh3/minimap2)
- [emtrey](https://github.com/rvolden/emtrey) ([go](https://golang.org/dl/))
- [blat source](https://users.soe.ucsc.edu/~kent/src/blatSrc35.zip) or [blat executable](http://hgdownload.soe.ucsc.edu/admin/exe/)
- [medaka](https://github.com/nanoporetech/medaka)

These will need to be put into your path.

Note: We plan to discontinue medaka use in the next version. Medaka is very slow, requires us to update the polishing models, and has only a marginal benefit on isoform accuracy. When we do this, we will discontinue the "PC" consensus mode

The following python libraries are required as well:

- [mappy](https://pypi.org/project/mappy/)
- [pyabPOA](https://pypi.org/project/pyabpoa/)

```bash
pip3 install mappy
pip3 install pyabpoa==1.0.5
```
should just work to install them. Make sure to get the right version (1.0.5) of pyabpoa

## Usage ##
```bash
python3 Mando.py [OPTIONS]
```

Running with default settings:
```bash
python3 Mando.py -p . -g gencodeV29.gtf -G hg38.fasta -f Consensus_reads.fofn
```
Everything else is optional:

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

  -m, --multi_exon_only
                        If used, Mandalorion will filter all single exon
                        isoforms

  -C CONSENSUSMODE, --consensusMode CONSENSUSMODE
                        Set to P or PC (deault = P). If P, only pyabpoa will
                        be used to make isoform consensus sequences. If PC,
                        medaka will be used as well. PC generates slightly
                        more accurate sequences for R2C2 data but is much slower.
                        Use P is analyzing PacBio data.

  -v, --version         Prints Mandalorion version
```

## Outputs ##

I consider the *Isoforms.filtered.fasta* file the main output of the Mandalorion pipeline. It contains the polished sequences of all isoforms Mandalorion considers very high confidence. Mandalorion also creates *Isoforms.filtered.clean.psl* and *Isoforms.filtered.clean.gtf* files which contain minimap2 alignments of those sequences that had small indels removed. These files are easy to upload to the UCSC Genome Browser to inspect. This version of Mandalorion now also generates a *Isoforms.filtered.clean.quant* file which contains the number of R2C2/PacBio reads that associate which each isoform for each given fasta files.

## Utils ##

These are the scripts used to do haplotype phasing and HLA analysis as well as preparing reads for Mandalorion and parsing output to fit the LRGASP consortium requirements. 
