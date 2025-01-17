import sys
import argparse
import numpy
import gzip

def argParser():
    parser = argparse.ArgumentParser(
        description='Groups isoforms into loci and matches them with genes if annotation gtf is given',
        add_help=True,
        prefix_chars='-',
    )

    parser.add_argument('-i', '--infile', type=str, action='store')
    parser.add_argument('-o', '--outfile', type=str, action='store')
    parser.add_argument('-g', '--genome_annotation', type=str, action='store')


    return vars(parser.parse_args())


args = argParser()
infile = args['infile']
outfile = args['outfile']
genome_annotation = args['genome_annotation']

out=open(outfile,'w')

def read_annotation(genome_annotation):
    geneDict={}
    geneDict['+']={}
    geneDict['-']={}
    coordDict={}
    coordDict['+']={}
    coordDict['-']={}
    if genome_annotation == 'None':
        return coordDict

    else:
        if genome_annotation.endswith('.gtf.gz'):
            print('\t\tgtf file ends on .gz and will be treated as gzipped')
            input=gzip.open(genome_annotation,'rt')
        elif genome_annotation.endswith('.gtf'):
            input=open(genome_annotation,'r')

        with input as f:
            for line in f:
                if line[0]!='#':
                    a=line.strip().split('\t')
                    chromosome=a[0]
                    element=a[2]
                    left=int(a[3])-1
                    right=int(a[4])
                    direction=a[6]
                    gene_name=a[8].split('gene_id "')[1].split('"')[0]
                    if element == 'exon':
                        if 'gene_name' in a[8]:
                            gene_name+='_'+a[8].split('gene_name "')[1].split('"')[0]
                        if gene_name not in geneDict[direction]:
                            geneDict[direction][gene_name]=[chromosome,[]]
                        geneDict[direction][gene_name][1].append((left,right))


        for direction in ['+','-']:
            print('\t\treading genes on',direction,'strand')
            total=len(geneDict[direction])
            current=0
            for gene,coords in geneDict[direction].items():
                current+=1
                print('\t\t'+str(current),'of',total,str(round((current/total)*100,2))+'%'+' '*40,end='\r')
                chromosome=coords[0]
                for exon in coords[1]:
                    start=exon[0]
                    end=exon[1]
                    if chromosome not in coordDict[direction]:
                        coordDict[direction][chromosome]={}
                    for i in range(start,end,2):
                        if i not in coordDict[direction][chromosome]:
                            coordDict[direction][chromosome][i]=set()
                        coordDict[direction][chromosome][i].add(gene)
            print('\n')
        return coordDict


def group_isoforms(infile,coordDict):
    outDict={}
    isoforms=[]
    new=False
    previous_chrom=''
    previous_start=0
    previous_end=0
    roots=set()
    locus=0
    covered=set()
    for direction in ['+','-']:
        print('\t\tgrouping isoforms  on',direction,'strand')
        outDict={}
        isoforms=[]
        new=False
        previous_chrom=''
        previous_start=0
        previous_end=0
        roots=set()
        locus=0
        isoform_count=0
        for line in open(infile):
            a=line.strip().split('\t')
            isodirection=a[8]
            if direction==isodirection:
                isoform_count+=1

        current=0
        for line in open(infile):
            a=line.strip().split('\t')
            isodirection=a[8]
            if direction==isodirection:
                current+=1
                print('\t\t'+str(current),'of',isoform_count,str(round((current/isoform_count)*100,2))+'%'+' '*40,end='\r')

                chrom=a[13]


                start=int(a[15])
                end=int(a[16])

                if chrom!=previous_chrom:
                    new=True
                else:
                    if start>previous_end:
                        new=True
                    else:
                        previous_end=max(end,previous_end)
                        isoforms.append(line)
                if new:
                    if isoforms:
                        locus = match_isoforms(isoforms,previous_chrom,previous_start,previous_end,direction,locus)
                    isoforms=[]
                    isoforms.append(line)
                    previous_chrom=chrom
                    previous_end=end
                    previous_start=start
                    new=False
        if isoforms:
            locus = match_isoforms(isoforms,previous_chrom,previous_start,previous_end,direction,locus)
        print('\n')

def match_isoforms(isoforms,previous_chrom,previous_start,previous_end,direction,locus):
    genes={}
    covered=set()
    for line in isoforms:
        a=line.strip().split('\t')
        blockstarts=numpy.array(a[20].split(',')[:-1],dtype=int)
        blockwidths=numpy.array(a[18].split(',')[:-1],dtype=int)
        for index in range(0,len(blockstarts),1):
            blockstart=blockstarts[index]
            blockend=blockstart+blockwidths[index]
            for i in range(blockstart,blockend,1):
                covered.add(i)


    for i in covered:
        if previous_chrom in coordDict[direction]:
            if i in coordDict[direction][previous_chrom]:
                for gene in coordDict[direction][previous_chrom][i]:
                    if gene not in genes:
                        genes[gene]=0
                    genes[gene]+=1

    geneList=[]
    geneSet=set()
    for gene,count in genes.items():
        geneList.append((count,gene))
        geneSet.add(gene)
    coords=[0]

    if geneList:
        best=sorted(geneList,reverse=True)[0][1]
    else:
        best=''
    locus+=1
    LocusName='Locus'+str(locus)
    if geneSet:
        GeneOverlaps=(',').join(geneSet)
    else:
        GeneOverlaps=''
    for line in isoforms:
        name=line.strip().split('\t')[9]
        out.write(name+'\t'+LocusName+'\t'+previous_chrom+'\t'+str(previous_start)+'\t'+str(previous_end)+'\t'+best+'\t'+GeneOverlaps+'\n')
    return locus


print('\t\tparsing genome annotation')
coordDict=read_annotation(genome_annotation)
print('\t\tprocessing isoforms')
group_isoforms(infile,coordDict)
