#!/usr/bin/env python3
# Christopher Vollmers
# Roger Volden

import sys
import os
import argparse
import mappy as mp
import time
from time import localtime, strftime


VERSION = 'v3.6.3 - I find this isoform vague and unconvincing'

parser = argparse.ArgumentParser(usage='\n\nRunning with default parameters:\n\npython3 Mando.py -p . -g gencodeV29.gtf -G hg38.fasta -f Consensus_reads_noAdapters_noPolyA_5->3.fofn\n')

parser.add_argument(
    '-p', '--path', type=str, help='Directory to put output files into'
)
parser.add_argument(
    '-u', '--upstream_buffer', type=str, default='10',
    help='''Defines upstream leniency window for
            polyA and TSS determination (default 10)'''
)
parser.add_argument(
    '-d', '--downstream_buffer', type=str, default='50',
    help='''Defines downstream leniency window for polyA
            and TSS determination (default 50)'''
)
parser.add_argument(
    '-g', '--genome_annotation', type=str, default='None',
    help='''Genome annotation file (gtf).
            Is used to identify individual annotated splice sites.
            If -W is set it is also used to white-list annotated polyA sites'''
)
parser.add_argument(
    '-G', '--genome_sequence', type=str, help='Genome file (fasta)'
)
parser.add_argument(
    '-r', '--minimum_ratio', type=str, default='0.01',
    help='''Proportion of reads that align to a locus
            required for an isoform (default 0.01)'''
)
parser.add_argument('-i', '--minimum_internal_ratio', type=str, default='1')
parser.add_argument(
    '-R', '--minimum_reads', type=str, default='3',
    help='''Minimum number of reads to make an isoform (default 3)'''
)

parser.add_argument(
    '-f', '--Consensus_reads', type=str,
    help='''Fasta/fastq file with R2C2/PacBio consensus reads,
            can be entered as a single file path,
            a comma separated list of  file paths,
            or a path to a file of filenames file (has to end on .fofn) that contains one file path per line'''
)

parser.add_argument(
    '-O', '--overhangs', type=str, default='0,40,0,40',
    help='''Defines bounds for unaligned bases on ends. Format:
            min5prime,max5prime,min3prime,max3prime (default 0,40,0,40)'''
)
parser.add_argument(
    '-t', '--minimap2_threads', type=str, default='4',
    help='Number of threads to use when running minimap and consensus calling (default 4)'
)

parser.add_argument(
    '-I', '--minimum_isoform_length', type=str, default='500',
    help='Minimum length in nt of isoforms that will be considered (default 500)'
)
parser.add_argument(
    '-n', '--minimum_feature_count', type=str, default='2',
    help='''Features (starts,ends, new splice sites) will be considered
            if they are present in this number of reads (default 2)'''
)
parser.add_argument(
    '-w', '--splice_site_window', type=str, default='1',
    help='''reads spliced within this number of nucleotides on each side
            of a splice site will be considered spliced at this site (default 1)'''
)
parser.add_argument(
    '-A', '--Acutoff', type=str, default='0.5',
    help='''Isoforms with A content of more than this cutoff in a 30nt
            window surrounding their polyA site will be discarded (default 0.5)'''
)
parser.add_argument(
    '-W', '--white_list_polyA', type=str, default='0',
    help='''If set, polyA sites that fall within +/-20nt of annotated transcript ends will not be filtered regardless of Acutoff set with -A.
            Annotated transcript ends will be  taken from annotation file given with -g.
            Only transcripts having one of the provided comma separated values in the line of their exon features will be used.
            E.g. setting -W to [SIRV,"basic"] will whitelist spike-in SIRV transcripts and transcripts considered "basic", i.e. high confidence full length, in the gencode annotation.
            Setting -W to [exon] should include all transcripts in the gtf file.
            This is feature only checks whether the provided values are in the line, not where they are.
            That means that setting -W to [chr1] will whitelist all transcripts on chromosome 1 but also transcripts with "chr1" in their name, so it's a bit dangerous'''
)

parser.add_argument(
    '-S', '--sam_file', type=str, default=False,
    help='''If given, Mandalorion will use this file instead of performing its own minimap2 alignment.
            Careful! If names don't line up between this and the fasta and fastq files everything breaks!'''
)

parser.add_argument(
    '-m', '--multi_exon_only', default='0',action='store_const', const='1',
    help='''If used, Mandalorion will filter all single exon isoforms'''
)


parser.add_argument(
    '-M', '--Modules', default='APSDCTFQ',
    help='''Defines what modules of Mandalorion will be run. By default this includes:
            A - Alignment,
            P - .sam to .clean.psl conversion,
            S - def splice sites,
            D - defining isoforms,
            C - Creating consensus sequences for those isoforms,
            F - Filtering isoforms,
            Q - Quantifying isoforms.
            Each module needs the output of the modules run before it to function properly.
            Running individual modules can be useful if you want to for example filter with different parameters without rerunning the whole pipeline. (default APSDCTFQ)'''
)


parser.add_argument(
    '-v', '--version', action='version', version=VERSION,
    help='Prints Mandalorion version'
)

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(0)
args = parser.parse_args()

path = args.path + '/'  # path where you want your output files to go
temp_path = path + '/tmp/'
upstream_buffer = args.upstream_buffer
downstream_buffer = args.downstream_buffer
genome_annotation = args.genome_annotation
genome_sequence = args.genome_sequence
minimum_ratio = args.minimum_ratio
minimum_internal_ratio = args.minimum_internal_ratio
minimum_reads = args.minimum_reads
fasta_files = args.Consensus_reads
overhangs = args.overhangs
minimap2_threads = args.minimap2_threads
minimum_isoform_length = args.minimum_isoform_length
window = args.splice_site_window
feature_count = args.minimum_feature_count
Acutoff = args.Acutoff
sam_file=args.sam_file
white_list_polyA=args.white_list_polyA
multi_exon_only=args.multi_exon_only
Modules=args.Modules

MandoPath = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/'


if '.fofn' in fasta_files:
    fastaList=[]
    for line in open(fasta_files):
        fasta=line.strip()
        fastaList.append(fasta)
else:
    fastaList=fasta_files.split(',')

if not os.path.isdir(path):
    os.system('mkdir %s' % path)
    command=open(path+'/Mando.log','w')
else:
    command=open(path+'/Mando.log','a')

command.write('\nMandalorion "'+ VERSION +'" was run on '\
               + strftime("%Y-%m-%d %H:%M:%S", localtime())+'\n'\
               + 'with the following parameters\n'\
               + str(args).replace('Namespace(','').replace(')','') + '\n')
command.close()
if not os.path.isdir(temp_path):
    os.system('mkdir %s' % temp_path)




minimap2 = 'minimap2'
emtrey = 'emtrey'

print('\n-------------------------------------------------------------\
       \nRunning Mandalorion - Isoform identification pipeline\
       \nVersion -', VERSION,'\
       \n-------------------------------------------------------------\n')


if not sam_file:
    sam_file = temp_path + '/mm2Alignments.sam'
    if 'A' in Modules:
        print('\n----------------------------\
               \nRunning Module A - Alignment\
               \n----------------------------\n')
        tempFasta=temp_path+'/Combined.fasta'
        out=open(tempFasta,'w')
        for fasta in fastaList:
            for name,seq,qual in mp.fastx_read(fasta):
                out.write('>%s\n%s\n' %(name,seq))

        out.close()
        os.system(
            '%s -G 400k --secondary=no -ax splice:hq --cs=long -uf -t %s %s %s > %s '
            % (minimap2, minimap2_threads, genome_sequence, tempFasta, sam_file)
        )

else:
    print('\n\n**********  sam file provided. Using those alignments instead of aligning the reads in fasta file(s)  **********\n\n')


psl_file = temp_path + '/mm2Alignments.psl'
clean_psl_file = temp_path + '/mm2Alignments.clean.psl'
if 'P' in Modules:
    print('\n----------------------------------------\
           \nRunning Module P - sam to Psl conversion\
           \n----------------------------------------\n')

    print('\tConverting sam output to psl format\n')
    os.system('%s -i %s > %s ' % (emtrey, sam_file, psl_file))
    print('\tCleaning psl file of small Indels\n')
    os.system('python3 %s/%s %s %s ' % (MandoPath,'clean_psl.py', psl_file, clean_psl_file))

if 'S' in Modules:
    print('\n----------------------------------------\
           \nRunning Module S - defining splice sites\
           \n----------------------------------------\n')
    start=time.time()
    os.system(
        'python3 %s/spliceSites.py -i %s -p %s -c %s -g %s -r %s -s %s -w %s -m %s -W %s'# -n %s'
        % (
            MandoPath,
            clean_psl_file,
            temp_path,
            '0.1',
            genome_annotation,
            'g',
            sam_file,
            window,
            feature_count,
            white_list_polyA#,
#            minimap2_threads
        )
    )
    print('duration',time.time()-start)
if 'D' in Modules:
    print('\n------------------------------------\
           \nRunning Module D - defining isoforms\
           \n------------------------------------\n')
# This script sort raw reads into isoform bins.
# The two number variables determine the window around TSS and TES
# in which read ends can fall and still be matched to the site.
    os.system(
        'python3 %s/defineAndQuantifyIsoforms.py -i %s -p %s -d %s -u %s -f %s -n %s -R %s'
        % (
            MandoPath,
            clean_psl_file,
            temp_path,
            downstream_buffer,
            upstream_buffer,
            fasta_files,
            feature_count,
            minimum_reads
        )
    )
if 'C' in Modules:
    print('\n-----------------------------------------------\
           \nRunning Module C - creating consensus sequences\
           \n-----------------------------------------------\n')
    os.system(
        'python3 %s/createConsensi.py -p %s -n %s'
        % (MandoPath,temp_path, minimap2_threads)
    )

    os.system('scp %s %s' % (temp_path + '/isoform_tmp.fasta',temp_path + 'Isoforms_full_length_consensus_reads.fasta'))

if 'F' in Modules:
    print('\n-------------------------------------\
           \nRunning Module F - filtering isoforms\
           \n-------------------------------------\n')
    os.system(
        'python3 %s/filterIsoforms.py \
            -p %s -i %s -r %s -R %s -n %s -G %s \
            -O %s -t %s -A %s -s %s -d %s -I %s -m %s -M %s 2> %s'
        % (
            MandoPath,
            temp_path,
            temp_path + '/Isoform_Consensi.fasta',
            minimum_ratio,
            minimum_reads,
            minimum_internal_ratio,
            genome_sequence,
            overhangs,
            minimap2_threads,
            Acutoff,
            window,
            downstream_buffer,
            minimum_isoform_length,
            MandoPath,
            multi_exon_only,
            temp_path + '/filter_reasons.txt',
        )
    )
    os.system('scp ' + temp_path + '/Isoforms.filtered.* ' + path)
if 'Q' in Modules:
    print('\n---------------------------------------\
           \nRunning Module Q - quantifying isoforms\
           \n---------------------------------------\n')
    os.system(
        'python3 %s/assignReadsToIsoforms.py -m %s -f %s'
        % (MandoPath,temp_path, fasta_files)
    )
    os.system('scp ' + temp_path + '/Isoforms.filtered.clean.quant ' + path)
