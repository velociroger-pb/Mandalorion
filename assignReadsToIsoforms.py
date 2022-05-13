import sys
import argparse
import mappy as mp


parser = argparse.ArgumentParser()
parser.add_argument('-m', '--mandalorion_output_folder', type=str)
parser.add_argument('-f', '--fasta_files', type=str, help='comma separate list of fasta file locations')



args = parser.parse_args()
mandalorion_folder=args.mandalorion_output_folder
fasta_files=args.fasta_files
filtered_isoforms=mandalorion_folder+'/Isoforms.filtered.clean.psl'
r2i=mandalorion_folder+'/reads2isoforms.txt'
outq=open(mandalorion_folder+'/Isoforms.filtered.clean.quant','w')

def read_fasta(inFile):
    '''Reads in FASTA files, returns a dict of header:sequence'''
    readDict = {}
    for name,seq,qual in mp.fastx_read(inFile):
         readDict[name] = ''

    return readDict

def read_filtered_isoforms(filtered_isoforms,r2i_dict,sampleList,readMapDict):
    for line in open(filtered_isoforms):
        a=line.strip().split('\t')
        isoform=a[9]
        quantDict={}
        for sample in sampleList:
            quantDict[sample]=0

        for name in r2i_dict[isoform]:
            sample=readMapDict[name]
            quantDict[sample]+=1

        outq.write(a[9]+'\t')
        for sample in sampleList:
            outq.write(str(quantDict[sample])+'\t')
        outq.write('\n')

def mapReadLocation(fastaList):
    sampleList=[]
    readMapDict={}
    for line in fastaList:
        location=line.strip()
        sampleList.append(location)
        reads=read_fasta(location)
        for name,seq,qual in mp.fastx_read(location):
            readMapDict[name]=location
    outq.write('Isoform\t')
    for sample in sampleList:
        outq.write(sample+'\t')
    outq.write('\n')
    return sampleList,readMapDict


def read_r2i(r2i):
    r2i_dict={}
    for line in open(r2i):
        a=line.strip().split('\t')
        read=a[0]
        isoform=a[1]
        if isoform not in r2i_dict:
            r2i_dict[isoform]=[]
        r2i_dict[isoform].append(read)
    return r2i_dict


if '.fofn' in fasta_files:
    fastaList=[]
    for line in open(fasta_files):
        fasta=line.strip()
        fastaList.append(fasta)
else:
    fastaList=fasta_files.split(',')





sampleList,readMapDict=mapReadLocation(fastaList)
r2i_dict=read_r2i(r2i)
read_filtered_isoforms(filtered_isoforms,r2i_dict,sampleList,readMapDict)
