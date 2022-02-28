# Mandalorion #
[![Github release](https://img.shields.io/github/tag/christopher-vollmers/Mandalorion-1.svg?label=Version)](https://github.com/christopher-vollmers/Mandalorion-1/tags)
[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](http://perso.crans.org/besson/LICENSE.html)


v3.6.3: I find this isoform vague and unconvincing.

Takes R2C2/C3POa and/or PacBio/ccs/lima data and defines high confidence isoform consensus sequences and alignments.
You can mix and match R2C2/PacBio reads and fasta/fastq files (quality scores are ignored).
Mandalorion is not tested for use with regular ONT reads. I recommend using FLAIR instead for those.

## PolyA tail removal ##

It is best to run Mandalorion with reads that have no adapter sequences (trimmed), no polyA tails, and in 5'->3' orientation. By default, C3POa (v2.3.0) does trimm reads but doesn't do polyA removal and, depending on the adapters used, may not put reads in 5'->3' orientation..   
Therefore the removePolyA.py script in utils/ should be used to get both PacBio (depending on preprocessing) and R2C2 reads ready.
The script first removes a fixed number of bases (can't be zero) from the ends of the read and then removes the polyA tail on the 3' end of the read. Depending on the polyA position it will also reorient the read in 5'->3' direction.

## Workflow

Mandalorion contains several modules that are run sequentially by default (APSDCFQ) to define, align, and quantify high-quality isoforms.

![Asset1](https://user-images.githubusercontent.com/28308271/156075738-51d545ca-0b4c-4ef5-b2ac-2b29f36bc552.png)

## Dependencies ##

- [minimap2](https://github.com/lh3/minimap2)
- [emtrey](https://github.com/rvolden/emtrey) ([go](https://golang.org/dl/))

These will need to be put into your path.

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


  -f CONSENSUS_READS, --Consensus_reads CONSENSUS_READS
                        Fasta/fastq file with R2C2/PacBio consensus reads, can be entered
                        as a single file path, a comma separated list of file
                        paths, or a path to a file of filenames file (has to
                        end on .fofn) that contains one file path per line.

  -O OVERHANGS, --overhangs OVERHANGS
                        Defines bounds for unaligned bases on ends. Format:
                        min5prime,max5prime,min3prime,max3prime (default
                        0,40,0,40)

  -t MINIMAP2_THREADS, --minimap2_threads MINIMAP2_THREADS
                        Number of threads to use when running minimap and
                        consensus calling (default 4)

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
                        default this includes:
                        A - Alignment,
                        P - .sam to .clean.psl conversion,
                        S - defining splice sites,
                        D - defining isoforms,
                        C - Creating consensus sequences for those isoforms
                        F - Filtering isoforms
                        Q - Quantifying isoforms.

                        Each module need the output of the modules run before it to
                        function properly. Running individual modules can be
                        useful if you want to for example filter with
                        different parameters without rerunning the whole
                        pipeline

  -m, --multi_exon_only
                        If used, Mandalorion will filter all single exon
                        isoforms

  -v, --version         Prints Mandalorion version
```

## Outputs ##

I consider the *Isoforms.filtered.fasta* file the main output of the Mandalorion pipeline. It contains the polished sequences of all isoforms Mandalorion considers very high confidence. Mandalorion also creates *Isoforms.filtered.clean.psl* and *Isoforms.filtered.clean.gtf* files which contain minimap2 alignments of those sequences that had small indels removed. These files are easy to upload to the UCSC Genome Browser to inspect. This version of Mandalorion now also generates a *Isoforms.filtered.clean.quant* file which contains the number of R2C2/PacBio reads that associate which each isoform for each given fasta files.

## Utils ##

These are the scripts used to do haplotype phasing and HLA analysis as well as preparing reads for Mandalorion and parsing output to fit the LRGASP consortium requirements.
