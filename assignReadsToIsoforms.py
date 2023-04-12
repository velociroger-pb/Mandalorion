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
outq=open(mandalorion_folder+'/Isoforms.filtered.clean.counts','w')
outtpm=open(mandalorion_folder+'/Isoforms.filtered.clean.tpm','w')

def read_fasta(inFile):
    '''Reads in FASTA files, returns a dict of header:sequence'''
    readDict = {}
    for name,seq,qual in mp.fastx_read(inFile):
         readDict[name] = ''

    return readDict

def read_filtered_isoforms(filtered_isoforms,r2i_dict,sampleList,readMapDict,isoformReadCounts,totalReadCounts):
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
        outtpm.write(a[9]+'\t')
        for sample in sampleList:
            value=quantDict[sample]
            isoformReads=isoformReadCounts[sample]
            totalReads=totalReadCounts[sample]
            outq.write(str(value)+'\t')
            outtpm.write(str(round((value/totalReads)*1000000,3))+'\t')
        outq.write('\n')
        outtpm.write('\n')


def mapReadLocation(fastaList):
    sampleList=[]
    readMapDict={}
    totalReadCounts={}
    for line in fastaList:
        location=line.strip()
        totalReadCounts[location]=0
        sampleList.append(location)
        reads=read_fasta(location)
        for name,seq,qual in mp.fastx_read(location):
            readMapDict[name]=location
            totalReadCounts[location]+=1
    outq.write('Isoform\t')
    outtpm.write('Isoform\t')
    for sample in sampleList:
        outq.write(sample+'\t')
        outtpm.write(sample+'\t')
    outq.write('\n')
    outtpm.write('\n')
    return sampleList,readMapDict,totalReadCounts


def read_r2i(r2i,readMapDict):
    r2i_dict={}
    isoformReadCounts={}
    for line in open(r2i):
        a=line.strip().split('\t')
        read=a[0]
        isoform=a[1]
        location=readMapDict[read]
        if location not in isoformReadCounts:
            isoformReadCounts[location]=0
        if isoform not in r2i_dict:
            r2i_dict[isoform]=[]
        r2i_dict[isoform].append(read)
        isoformReadCounts[location]+=1
    return r2i_dict,isoformReadCounts


if '.fofn' in fasta_files:
    fastaList=[]
    for line in open(fasta_files):
        fasta=line.strip()
        fastaList.append(fasta)
else:
    fastaList=fasta_files.split(',')




print('\t\tmap reads to files')
sampleList,readMapDict,totalReadCounts=mapReadLocation(fastaList)
print('\t\tread to isoform info')
r2i_dict,isoformReadCounts=read_r2i(r2i,readMapDict)
print('\t\tquantify isoforms')
read_filtered_isoforms(filtered_isoforms,r2i_dict,sampleList,readMapDict,isoformReadCounts,totalReadCounts)
