import sys
import argparse
import mappy as mp


parser = argparse.ArgumentParser()
parser.add_argument('-m', '--mandalorion_output_folder', type=str)
parser.add_argument('-f', '--fasta_files', type=str, help='comma separate list of fasta file locations')



args = parser.parse_args()
mandalorion_folder=args.mandalorion_output_folder
fasta_files=args.fasta_files
filtered_isoforms=mandalorion_folder+'/Isoform_Consensi_filtered.aligned.out.clean.psl'
isoform_long_names=mandalorion_folder+'/Isoform_long_names.txt'
isoform_list=mandalorion_folder+'/isoform_list'
read_location=mandalorion_folder+'/read_locations.txt'
out=open(mandalorion_folder+'/reads2isoforms.txt','w')
outq=open(mandalorion_folder+'/quantifiedIsoforms.txt','w')

def read_gtf_file(gtf_file):
    gene_dict={}
    for line in open(gtf_file):
        a=line.strip().split('\t')
        if len(a)>6:

            type1=a[2]
            info=a[8]
            if type1=='gene':
                gene_id=info.split('gene_id "')[1].split('"')[0]
                if 'SIRV' in gene_id:
                    gene_symbol==gene_id
                else:
                    gene_symbol=info.split('gene_name "')[1].split('"')[0]
                gene_dict[gene_symbol]=(gene_id,a[3],a[4])
    return gene_dict


def read_sqanti_classification(sqanti_file):
    gene_dict={}
    for line in open(sqanti_file):
        a=line.strip().split('\t')
        isoform=a[0]
        gene=a[6]
        chromosome=a[1]
        gene_dict[isoform]=(gene,chromosome)
    return gene_dict


def read_fasta(inFile):
    '''Reads in FASTA files, returns a dict of header:sequence'''
    readDict = {}
    for name,seq,qual in mp.fastx_read(inFile):
         readDict[name] = ''

    return readDict


def read_isoform_long_names(isoform_long_names):
    short2longDict={}
    for line in open(isoform_long_names):
        a=line.strip().split('\t')
        long_name=('_').join(a[0].split('_')[:-1])
        short_name=a[1]
        short2longDict[short_name]=long_name
    return short2longDict

def read_filtered_isoforms(filtered_isoforms,short2longDict,long2locationDict,sampleList,readMapDict):#,geneDict,geneSymbols):
    for line in open(filtered_isoforms):
        a=line.strip().split('\t')
        short_name=('_').join(a[9].split('_')[:-1])
        long_name=short2longDict[short_name]
        location=long2locationDict[long_name]
        quantDict={}
#        gene,chromosome=geneDict[a[9]]
#        if gene in geneSymbols:
#            gene_symbol,start,end=geneSymbols[gene]
#        else:
#            gene_symbol,start,end='-','-','-'
#
        for sample in sampleList:
            quantDict[sample]=0
        for name,seq,qual in mp.fastx_read(location):
            out.write(name+'\t'+a[9]+'\n')
            sample=readMapDict[name]
            quantDict[sample]+=1

        outq.write(a[9]+'\t')
        for sample in sampleList:
            outq.write(str(quantDict[sample])+'\t')
        outq.write('\n')

def read_isoform_list(isoform_list):
    long2locationDict={}
    for line in open(isoform_list):
        a=line.strip().split('\t')
        location=a[0]
        long_name=a[2]
        long2locationDict[long_name]=location
    return long2locationDict


def mapReadLocation(fastaList):
    sampleList=[]
    readMapDict={}
    for line in fastaList:
        print(line)
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


if '.fofn' in fasta_files:
    fastaList=[]
    for line in open(fasta_files):
        fasta=line.strip()
        fastaList.append(fasta)
else:
    fastaList=fasta_files.split(',')

sampleList,readMapDict=mapReadLocation(fastaList)
print(sampleList)
short2longDict=read_isoform_long_names(isoform_long_names)
long2locationDict=read_isoform_list(isoform_list)
read_filtered_isoforms(filtered_isoforms,short2longDict,long2locationDict,sampleList,readMapDict)#,geneDict,geneSymbols)
