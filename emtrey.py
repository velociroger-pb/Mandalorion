import argparse
import multiprocessing as mp
import gc
import re
import numpy as np
import mappy


parser = argparse.ArgumentParser()

parser.add_argument('--threads','-t',default='50',type=str,action='store',help='number of threads to be used')
parser.add_argument('--outputFile','-o',type=str,action='store',help='output file goes here')
parser.add_argument('--inputFile','-i',type=str,action='store',help='path to input file')
parser.add_argument('--batch','-b',default='100000',type=str,action='store',help='lines are processed in batches of this size')
parser.add_argument('--mando','-m',default=False,action='store_true',help='generates invalid psl file for internal Mandalorion use')



args = parser.parse_args()

threads=int(args.threads)
inFile=args.inputFile
outputFile=args.outputFile
mando=args.mando
batch=int(args.batch)

out=open(outputFile,'w')

mm=True

def parseLine(a,qSize):
    name = a[0]
    tstart=int(a[3])-1
    bitwise=format(int(a[1]), "b")[::-1]
    if len(bitwise)>4:
        direction=bitwise[4]
    else:
        direction='0'

    if direction=='0':
        strand='+'
    else:
        strand='-'

    cstr=a[5]
    sequence=a[9]
    if strand=='-':
        sequence=mappy.revcomp(sequence)

    p = re.compile(r'([MIDNSHP=X])')
    splitCstr = [i+j for i,j in zip(p.split(cstr)[::2], p.split(cstr)[1::2])]
    blocksizes,qstarts,tstarts=[],[],[tstart]
    qstart,M,I,nI,D,nD,N,S,H,EQ,X,qend=0,0,0,0,0,0,0,0,0,0,0,0
    for i in range(len(splitCstr)):
        number=int(splitCstr[i][:-1])
        letter=splitCstr[i][-1]

        if letter in 'SH':
            if i==0:
                qstart=number
            elif i==len(splitCstr)-1:
                qend=number
        if i == 0:
            qstarts.append(qstart)
        if letter == 'M':
            m=number
            M+=m
            blocksizes.append(m)
            qstarts.append(m+qstarts[-1])
            tstarts.append(m+tstarts[-1])
        if letter == 'I':
            i=number
            I+=i
            nI+=1
            qstarts[-1]+=i
        if letter == 'D':
            d=number
            D+=d
            nD+=1
            tstarts[-1]+=d
        if letter == 'N':
            n=number
            N+=n
            tstarts[-1]+=n
        if letter =='S':
            s=number
            S+=s
        if letter =='H':
            h=number
            H+=h
        if letter =='=':
            eq=number
            EQ+=eq
        if letter =='X':
            x=number
            X+=x



    ID = I+D
    sLen = M + I + S + H + EQ + X
    consumeRef = M + D + N + EQ + X
    tend = tstart + consumeRef
    if qend == 0:
        end = sLen
    else:
        end = sLen - qend

    if len(qstarts) > 0:
        qstarts = qstarts[:-1]

    if len(tstarts) > 0:
        tstarts = tstarts[:-1]

    blockCount = len(blocksizes)

    NM, ambig, matches, mismatch = 0,0,0,0
    if mm:
        for col in a[9:]:
            if 'NM:i:' in col:
                NM = int(col.split(':')[2])
            if 'nn:i:' in col:
                ambig = int(col.split(':')[2])
            if 'ts:A:' in col:
                newStrand = col.split(':')[2]
                if newStrand == '-' and strand == '+':
                    strand = "-"
                elif newStrand == "-" and strand == "-":
                    strand = "+"
            if 'cs:Z:' in col:
                cs = col.split(':')[2]
        mismatch = NM - ID - ambig
        if mismatch < 0:
            mismatch = 0

        matches = M - mismatch
        accuracy = matches/(matches+mismatch+ID+ambig)
    elif EQ != 0:
        matches, mismatch = EQ, X
    else:
        mismatch, matches = 0, M

    bSize=(',').join(np.array(blocksizes,dtype=str))+','
    qSt=(',').join(np.array(qstarts,dtype=str))+','
    tSt=(',').join(np.array(tstarts,dtype=str))+','
    pslLine=f'{matches}\t{mismatch}\t0\t{N}\t{nI}\t{I}\t{nD}\t{D}\t{strand}\t{name}\t{sLen}\t{qstart}\t{end}\t{a[2]}\t{qSize}\t{tstart}\t{tend}\t{blockCount}\t{bSize}\t{qSt}\t{tSt}'
    if mando:
        pslLine+=f'\t{accuracy}\t{cs}\t{sequence}\n'
    else:
        pslLine+='\n'

    return pslLine

def processSamBatch(reads,chromosomes):
    results={}
    pool = mp.Pool(processes=threads)
    counter=0
    for read in reads:
        counter+=1
        a=read.strip().split('\t')
        if a[2]=='*':
            continue
        else:
            chrLength=chromosomes[a[2]]
            results[counter]=pool.apply_async(parseLine,[a,chrLength])

    pool.close()
    pool.join()
    gc.collect()
    for counter in results:
        newLine=results[counter].get()
        out.write(newLine)

def delegatingPslConversion(samFile,batch):
    target=batch
    chromosomes={}
    counter=0
    reads=[]
    with open(samFile,'r') as f:
        for line in f:
            if not line.startswith('@'):
                counter+=1
                if counter%target==0:
                    print('\t\tprocessing read', counter,' '*20,end='\r')
                    processSamBatch(reads,chromosomes)
                    reads=[]

                reads.append(line.rstrip())
            else:
                if line.startswith('@SQ'):
                    a=line.rstrip().split('\t')
                    chr=a[1].split(':')[1]
                    length=int(a[2].split(':')[1])
                    chromosomes[chr]=length
    processSamBatch(reads,chromosomes)


def main():
    delegatingPslConversion(inFile,batch)
    out.close()

main()
