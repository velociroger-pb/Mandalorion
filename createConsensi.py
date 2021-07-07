#!/usr/bin/env python3
# Christopher Vollmers
# Roger Volden

import sys
import os
import numpy as np
import argparse
import multiprocessing as mp
import time
import signal
import mappy as mm
import pyabpoa as poa
poa_aligner = poa.msa_aligner()


def revComp(seq):
    '''Returns the reverse complement of a seq'''
    bases = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', '-': '-'}
    return ''.join([bases[x] for x in list(seq)])[::-1]


def argParser():
    '''Parses command line arguments'''
    parser = argparse.ArgumentParser(
        description='Makes consensus sequences from R2C2 reads.',
        add_help=True, prefix_chars='-')
    parser.add_argument(
        '--path', '-p', type=str, action='store', default=os.getcwd(),
        help='Directory where all the files are/where they will end up.\
              Defaults to your current directory.')
    parser.add_argument('--subsample', '-s', type=int, action='store')
    parser.add_argument('--numThreads', '-n', type=int, action='store')
    parser.add_argument('--consensusMode', '-C', type=str, default='P', action='store')
    parser.add_argument(
        '--config', '-c', type=str, action='store', default='',
        help='If you want to use a config file to specify paths to\
              programs, specify them here. Use for poa, racon, water,\
              blat, and minimap2 if they are not in your path.')
    return vars(parser.parse_args())


def configReader(configIn):
    '''Parses the config file.'''
    progs = {}
    for line in open(configIn):
        if line.startswith('#') or not line.rstrip().split():
            continue
        line = line.rstrip().split('\t')
        progs[line[0]] = line[1]
    # should have minimap, racon, consensus, blat, and emtrey
    possible = set(['racon', 'medaka'])
    inConfig = set()
    for key in progs.keys():
        inConfig.add(key)
    # check for missing programs
    # if missing, default to path
    for missing in possible - inConfig:
        if missing == 'consensus':
            path = 'consensus.py'
        else:
            path = missing
        progs[missing] = path
        sys.stderr.write('Using ' + str(missing)
                         + ' from your path, not the config file.\n')
    return progs


args = argParser()
path = args['path']
temp_folder = path + '/mp'
subsample = args['subsample']
numThreads = args['numThreads']
consensusMode = args['consensusMode']

if args['config']:
    progs = configReader(args['config'])
    racon = progs['racon']
    medaka = progs['medaka']
else:
    racon = 'racon'
    medaka = 'medaka_consensus'

def simplify(infile, outfile, namefile):
    isoforms = read_fasta(infile)
    out = open(outfile, 'w')
    out1 = open(namefile, 'w')
    counter = 0
    for isoform, sequence in isoforms.items():
        if len(sequence) > 0:
            counter += 1
            out.write('>Isoform' + '_' + str(counter) + '_' + isoform.split('_')[-1] + '\n' + sequence + '\n')
            out1.write(isoform + '\tIsoform' + '_' + str(counter) + '\n')
    out.close()
    out1.close()


def read_fastq_file(seq_file):
    '''
    Takes a FASTQ file and returns a list of tuples
    In each tuple:
        name : str, read ID
        seed : int, first occurrence of the splint
        seq : str, sequence
        qual : str, quality line
        average_quals : float, average quality of that line
        seq_length : int, length of the sequence
    '''
    read_list1 = []
    read_list2 = []

    burn = False
    previous = ''
    for fullName, seq, qual in mm.fastx_read(seq_file):
        name_root = fullName.split('_')
        name = name_root[0]
        number = name_root[1]

        if previous != name:
            burn = False
            previous = name

        if number == '0':
            burn = True

        seq_length = len(seq)

        if burn or 'I' in number:
            read_list2.append((name, seq, qual, seq_length))
        else:
            read_list1.append((name, seq, qual, seq_length))

    return read_list1, read_list2


def read_fasta(inFile):
    '''Reads in FASTA files, returns a dict of header:sequence'''
    reads = {}
    for read in mm.fastx_read(inFile, read_comment=False):
        reads[read[0]] = read[1]
    return reads



def determine_consensus(name, fasta, fastq, counter):
    '''Aligns and returns the consensus'''
    corrected_consensus = ''
    repeats = '0'
    fasta_reads = []
    fasta_read_dict = read_fasta(fasta)
    for read, seq in fasta_read_dict.items():
        fasta_reads.append((read, seq))
    repeats = str(len(fasta_reads))

    poa_cons = temp_folder + '/consensus.'+counter+'.fasta'
    if 'P' in consensusMode:
        indeces = np.random.choice(np.arange(0, len(fasta_reads)),
                                   min(len(fasta_reads), 100), replace=False)
        subsample_fasta_reads = []
        for index in indeces:
            subsample_fasta_reads.append(fasta_reads[index])

        first = subsample_fasta_reads[0][1]
        sequences=[]
        mm_align = mm.Aligner(seq=first, preset='map-ont')
        for read,sequence in subsample_fasta_reads:
            for hit in mm_align.map(sequence):
                 if hit.is_primary:
                     if hit.strand==1:
                         sequences.append(sequence)
                     elif hit.strand==-1:
                         sequences.append(mm.revcomp(sequence))


        res = poa_aligner.msa(sequences, out_cons=True, out_msa=False)
        if len(sequences)<=2:
            consensus_sequence = sequences[0]
        elif not res.cons_seq:
            consensus_sequence = sequences[0]
        else:
            consensus_sequence = res.cons_seq[0]
    else:
        consensus_sequence=fasta_reads[0][1]


    out_cons_file = open(poa_cons, 'w')
    out_cons_file.write('>Consensus\n' + consensus_sequence + '\n')
    out_cons_file.close()
    final=poa_cons
    corrected_consensus = consensus_sequence


    if 'C' in consensusMode:
        fastq_reads_full, fastq_reads_partial = read_fastq_file(fastq)
        fastq_reads = fastq_reads_full + fastq_reads_partial
        if len(fastq_reads) > 0:
            out_Fq = temp_folder + '/' + counter + '_subsampled.fastq'
            out = open(out_Fq, 'w')

            if len(fastq_reads_full) < subsample:
                subsample_fastq_reads = fastq_reads
            else:
                indeces = np.random.choice(
                    np.arange(0, len(fastq_reads_full)),
                    min(len(fastq_reads_full), subsample), replace=False)
                subsample_fastq_reads = []
                for index in indeces:
                    subsample_fastq_reads.append(fastq_reads_full[index])

            subread_counter = 0

            subsample_fastq_reads_numbered=[]
            for read in subsample_fastq_reads:
                subread_counter += 1
                out.write('@' + read[0] + '_' + str(subread_counter) + '\n'
                      + read[1] + '\n+\n' + read[2] + '\n')
                subsample_fastq_reads_numbered.append((read[0] + '_' + str(subread_counter), read[1], read[2], read[3]))
            out.close()
#            subsample_fastq_reads=list(subsample_fastq_reads_numbered)

#            overlap = temp_folder +'/overlaps.'+counter+'.paf'
#            overlap_fh=open(overlap,'w')
#            mm_align = mm.Aligner(seq=consensus_sequence, preset='map-ont')
#            for fqName,sequence,q,le in subsample_fastq_reads:
#                for hit in mm_align.map(sequence):
#                    if hit.is_primary:
#                        overlap_fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
#                            fqName, str(len(sequence)), hit.q_st, hit.q_en,
#                            hit.strand, 'Consensus', hit.ctg_len, hit.r_st,
#                            hit.r_en, hit.mlen, hit.blen, hit.mapq))

#            overlap_fh.close()

            output_cons = temp_folder + '/corrected_consensus.'+counter+'.fasta'
#            os.system('%s -q 5 -t 1 --no-trimming %s %s %s >%s 2>./racon_messages.txt' \
#                       %(racon,out_Fq, overlap, poa_cons, output_cons))
#            final=output_cons


            reads = read_fasta(final)
            for read in reads:
                corrected_consensus = reads[read] # if no read in file, corrected_consensus from pyabpoa output is used implicitly


            forMedaka = open(output_cons,'w')
            forMedaka.write('>Corrected_Consensus\n'+corrected_consensus+'\n')
            forMedaka.close()

            os.system('mkdir ' + temp_folder + '/' + counter)
            os.system('%s -f -i %s -d %s -o %s > %s_medaka_messages.txt 2>&1'
                      % (medaka, out_Fq, final,
                         temp_folder + '/' + counter, temp_folder + '/' + counter))
            final = temp_folder + '/' + counter + '/consensus.fasta'
            reads = read_fasta(final)
            for read in reads:
                corrected_consensus = reads[read]  # if no read in file, corrected_consensus from racon output is used implicitly

    combined_consensus_file = open(path + '/mp/' + counter + '.fasta', 'w')
    combined_consensus_file.write('>' + name + '_' + repeats + '\n'
                                  + corrected_consensus + '\n')
    combined_consensus_file.close()

def process_batch(batch):
     original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
     pool = mp.Pool(processes=numThreads)
     signal.signal(signal.SIGINT, original_sigint_handler)
     for name,fasta,fastq,sub_counter in batch:
         results=pool.apply_async(determine_consensus, [name, fasta, fastq, sub_counter])#.get()
     pool.close()
     pool.join()

def main():
    print('\tremoving files from previous run')
    os.system('rm -r ' + path + '/mp')
    print('\tdone removing files')
    os.system('mkdir ' + path + '/mp')
    if os.path.exists(path + '/Isoform_Consensi.fasta'):
        os.system('rm ' + path + '/Isoform_Consensi.fasta')
    os.system('touch {0}'.format(path + '/Isoform_Consensi.fasta'))
    batchSize=int(numThreads*20)
    counter = 0
    counter_list = []
    a=True
    if a:
#    try:
        generate=0
        for line in open(path + '/isoform_list'):
            generate+=1
        print('\tmaking consensus sequences of', generate, 'isoforms using consensusMode '+consensusMode)

        batch=[]
        for line in open(path + '/isoform_list'):
#          if 'SIRV4' in line:
            counter += 1
            fasta = line.split('\t')[0]
            fastq = line.split('\t')[1]
            name = line.split('\t')[2].strip()
            batch.append((name,fasta,fastq,str(counter)))
            counter_list.append(str(counter))
            if len(batch)==batchSize:
                process_batch(batch)
                batch=[]
                print('\tfinished consensus sequences of', len(counter_list), 'isoforms', end='\r')
        process_batch(batch)
        print('\tfinished consensus sequences of all isoforms')

#    except KeyboardInterrupt:
#        print('Caught KeyboardInterrupt, terminating workers')
#        pool.terminate()
#        pool.join()

    print('\tcombining temp files')
    combined_consensus_file = open(path + '/Isoform_Consensi.fasta', 'w')
    for counter in counter_list:
        if os.path.exists(path + '/mp/' + counter + '.fasta'):
            for lines in open(path + '/mp/' + counter + '.fasta'):
                combined_consensus_file.write(lines)
    combined_consensus_file.close()
    temp_fasta = path + '/isoform_tmp.fasta'
    simplify(path + '/Isoform_Consensi.fasta', temp_fasta, path + '/Isoform_long_names.txt')

main()
