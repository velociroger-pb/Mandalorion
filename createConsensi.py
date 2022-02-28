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
    parser.add_argument('--numThreads', '-n', type=int, action='store')

    return vars(parser.parse_args())




args = argParser()
path = args['path']
temp_folder = path + '/mp'
numThreads = args['numThreads']

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

def read_fasta(inFile):
    '''Reads in FASTA files, returns a dict of header:sequence'''
    reads = {}
    for read in mm.fastx_read(inFile, read_comment=False):
        reads[read[0]] = read[1]
    return reads



def determine_consensus(name, fasta, counter):
    '''Aligns and returns the consensus'''
    corrected_consensus = ''
    repeats = '0'
    fasta_reads = []
    fasta_read_dict = read_fasta(fasta)
    for read, seq in fasta_read_dict.items():
        fasta_reads.append((read, seq))
    repeats = str(len(fasta_reads))

    poa_cons = temp_folder + '/consensus.'+counter+'.fasta'
    indeces = np.random.choice(np.arange(0, len(fasta_reads)),
                                   min(len(fasta_reads), 100), replace=False)
    subsample_fasta_reads = []
    for index in indeces:
        subsample_fasta_reads.append(fasta_reads[index])
    first = subsample_fasta_reads[0][1]
    sequences=[]
    seq_lengths=[]
    mm_align = mm.Aligner(seq=first, preset='map-ont')
    for read,sequence in subsample_fasta_reads:
        seq_lengths.append(len(sequence))
        for hit in mm_align.map(sequence):
             if hit.is_primary:
                 if hit.strand==1:
                     sequences.append(sequence)
                 elif hit.strand==-1:
                     sequences.append(mm.revcomp(sequence))

    if np.median(seq_lengths)>14000:
        print('long')
        consensus_sequence = sequences[0]
    else:
        res = poa_aligner.msa(sequences, out_cons=True, out_msa=False)
        if len(sequences)<=2:
            consensus_sequence = sequences[0]
        elif not res.cons_seq:
            consensus_sequence = sequences[0]
        else:
            consensus_sequence = res.cons_seq[0]

    out_cons_file = open(poa_cons, 'w')
    out_cons_file.write('>Consensus\n' + consensus_sequence + '\n')
    out_cons_file.close()
    final=poa_cons
    corrected_consensus = consensus_sequence

    combined_consensus_file = open(path + '/mp/' + counter + '.fasta', 'w')
    combined_consensus_file.write('>' + name + '_' + repeats + '\n'
                                  + corrected_consensus + '\n')
    combined_consensus_file.close()

def process_batch(batch):
     original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
     pool = mp.Pool(processes=numThreads)
     signal.signal(signal.SIGINT, original_sigint_handler)
     for name,fasta,sub_counter in batch:
         results=pool.apply_async(determine_consensus, [name, fasta, sub_counter])#.get()
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
        print('\tmaking consensus sequences of', generate, 'isoforms')

        batch=[]
        for line in open(path + '/isoform_list'):
#          if 'SIRV4' in line:
            counter += 1
            fasta = line.split('\t')[0]
            name = line.split('\t')[1].strip()
            batch.append((name,fasta,str(counter)))
            counter_list.append(str(counter))
            if len(batch)==batchSize:
                process_batch(batch)
                batch=[]
                print('\tfinished consensus sequences of', len(counter_list), 'isoforms', end='\r')
        process_batch(batch)
        print('\tfinished consensus sequences of all isoforms                        ')

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
