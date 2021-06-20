import mappy
import os
import sys
import argparse


parser = argparse.ArgumentParser()

parser.add_argument('--infile_ccs', '-i', type=str, action='store', help='ccs reads in fastq or fastq.gz format')
parser.add_argument('--infile_subreads', '-s', type=str, action='store', help='takes fastq or fastq.gz files. bam files have to be converted first. Multiple files can be given as comma separated list. Subreads will be subsampled to 10 subreads per ccs read.')
parser.add_argument('--outfile_root', '-o', type=float, action='store', help='reads and subread output files will use this root')

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(0)
args = parser.parse_args()
infile_ccs = args.infile_ccs
infile_subreads = args.infile_subreads
outfile_root = args.outfile_root


outFile_fasta=open(outFile_root+'.fasta','w')
outFile_subread_fastq=open(outFile_root+'.subreads.fastq','w')

def qual2numbers(letter):
    number=ord(letter)-33
    return number


def readFastq(inFile_subreads,outFile_subreads_fastq):
    runningNumber=1
    coverageDict={}
    previousName=''
    previousLength=0
    rootName=''
    for subreadFile in inFile_subreads.split(','):
        for name,seq,qual in mappy.fastx_read(subreadFile):
            length=len(seq)

            splitName=name.split('/')

            rootName='PB-'+splitName[0].replace('_','-')+'-'+splitName[1]
            if rootName==previousName:
                runningNumber+=1
            else:
                coverageDict[previousName]=(runningNumber,previousLength)
                runningNumber=1
                previousLength=0

            previousName=rootName
            previousLength+=length

            newName=rootName+'_'+str(runningNumber)
            if runningNumber<=10:
                outFile_subread_fastq.write('@'+newName+'\n'+seq+'\n+\n'+qual+'\n')

    coverageDict[previousName]=(runningNumber,previousLength)
    return coverageDict


def readFasta(inFile_ccs,coverageDict,outFile_fasta):
    for name,seq,qual in mappy.fastx_read(inFile_ccs):
        length=len(seq)
        avgQ=15
        splitName=name.split('/')
        rootName='PB-'+splitName[0].replace('_','-')+'-'+splitName[1]
        coverage,rawLength=coverageDict[rootName]
        newName=rootName+'_'+str(avgQ)+'_'+str(rawLength)+'_'+str(coverage)+'_'+str(length)+'_'+str(length)
        outFile_fasta.write('>'+newName+'\n'+seq+'\n')


def main():
    coverageDict=readFastq(inFile_subreads,outFile_subread_fastq)
    readFasta(inFile_ccs,coverageDict,outFile_fasta)

main()


