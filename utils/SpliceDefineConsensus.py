#!/usr/bin/env python3
# Christopher Vollmers
# Roger Volden

import sys
import argparse
import re
import os
import numpy as np
import mappy as mp
import pyabpoa as poa
poa_aligner = poa.msa_aligner()


def get_parsed_files(out_tmp,chrom_list):
    roots=set()
    for file1 in os.listdir(out_tmp):
        if os.path.isfile(out_tmp+'/'+file1):
            if '.fasta' in file1:
                root=file1.split('.fasta')[0]
                chrom,start,end = root.split('~')
                chrom_list.add(chrom)
                roots.add(root)
    return chrom_list,roots

def readSAM(inFile,target_chrom):
    '''
    Reads the sam file, parses MD and CIGAR strings, and prints the read
    name, percent identity, number of matches, number of mismatches,
    and the number of indels.
    '''
    used = set()
    identities=[]
    lengths=[]
    readLengths=[]
    not_aligned=[]
    readAccuracy={}
    CigarDict={}

    for line in open(inFile):
        if line.startswith('@'):
            continue
        line = line.rstrip().split('\t')
        chromosome=line[2]
        if chromosome==target_chrom:
            name = line[0]
            readLength=len(line[9])
            chromosome=line[2]
            if line[1] == '4': # unaligned
                not_aligned.append(readLength)
            if line[1] == '4' or name in used: # unaligned
                continue
            readDirection=line[1]
            coverage=1
            if coverage:
                used.add(name)
                for column in line:
                    if column.startswith('NM:i:'):
                        NM = int(column[5:])
                    if column.startswith('nn:i:'):
                        ambig = int(column[5:])
                    if column.startswith('cs:Z:'):
                        CS=column[5:]
                CIGAR=CS+','+line[3]+','+line[5]+','+name
                denominator, M, indel = parseCIGAR(line[5])
                CigarDict[name]=CIGAR
                mismatch = NM - indel - ambig
                matches = M - mismatch
                identity = matches / denominator

                Length=round(M,-2)
                if name not in readAccuracy:
                    readAccuracy[name]=identity
    return readAccuracy,CigarDict

def parseCIGAR(cstr):
    '''
    Goes through the CIGAR string by pairing numbers with the letters.
    Counts the number of matches/mismatches and indels.
    Returns the total number of all, matches/mismatches (M), and indels.
    '''
    p = re.compile(r'([MIDNSHP=X])')
    splitCstr = [i+j for i,j in zip(p.split(cstr)[::2], p.split(cstr)[1::2])]
    total, M, indel = 0, 0, 0
    for i in range(len(splitCstr)):
        if splitCstr[i][-1] in 'MID':
            total += int(splitCstr[i][:-1])
        if splitCstr[i][-1] == 'M':
            M += int(splitCstr[i][:-1])
        if splitCstr[i][-1] in 'ID':
            indel += int(splitCstr[i][:-1])
    return total, M, indel

def getCSaroundSS(string,start,end):
    '''
    Goes through the CIGAR string by pairing numbers with the letters.
    Counts the number of matches/mismatches and indels.
    Returns the total number of all, matches/mismatches (M), and indels.
    '''
    CS,genomePosition, cstr, name = string.split(',')
    genomePosition=int(genomePosition)-1

    status=''
    p = re.compile(r'([\\=\\+\-\\*\\~])')
    splitCstr = [i+j for i,j in zip(p.split(CS)[1::2], p.split(CS)[2::2])]
    record=[]
    spliceIndex=False
    for i in range(len(splitCstr)):
        sub=splitCstr[i]
        status=sub[0]
        entry=sub[1:]
        if status == '=':
            for item in entry:
                genomePosition+=1
                record.append((status,item))
                if start<=genomePosition<=end:
                    spliceIndex=len(record)
        if status == '+':
            for item in entry:
                record.append((status,item))
        if status == '-':
            for item in entry:
                genomePosition+=1
                record.append((status,item))
                if start<=genomePosition<=end:
                    spliceIndex=len(record)
        if status == '*':
            for item in entry[::2]:
                genomePosition+=1
                record.append((status,item))
                if start<=genomePosition<=end:
                    spliceIndex=len(record)
        if status == '~':
            intronLength=int(entry[2:-2])
            genomePosition+=intronLength
            record.append('|'+entry+'|')
            if start<=genomePosition<=end:
                  spliceIndex=len(record)
    bases='nnnn'
    left=''
    right=''

    if spliceIndex:
        for index in range(max(spliceIndex-10,0),min(spliceIndex+10,len(record))):
            item=record[index]
            if '|' in item:
               bases=item[1:3]+item[-3:-1]
               left = record[index-5:index]
               right = record[index+1:index+6]
    return bases,left,right

def scan_for_best_bin(entry, dist_range, density_dict, peak_areas, chrom, side):
    best_extra_list, peak_center = 0, 0
    cov_area = {}
    best_dirs = {'+': 0, '-': 0}
    best_names = []
    for x in dist_range:
        extra_list_exp = 0
        dirs = {'+': 0, '-': 0}
        names=[]
        cov_set = {}
        called = False
        for y in dist_range:
            if entry + x + y in peak_areas[chrom][side]:
                called = True
        if called:
            continue
        for y in dist_range:
            if entry + x + y in density_dict:
                for item in density_dict[entry + x + y]:
                    extra_list_exp += 1
                    dirs[item[4]] += 1
                    names.append(item[0])
                    for covered_position in item[3]:
                        if covered_position not in cov_set:
                            cov_set[covered_position] = 0
                        cov_set[covered_position] += 1

        if extra_list_exp > best_extra_list:
            best_extra_list = extra_list_exp
            best_names=names
            peak_center = entry + x
            cov_area = cov_set
            best_dirs = dirs

    return best_extra_list, peak_center, cov_area, best_dirs,best_names


def determine_cov(cov_area, chrom, reverse, peak_center, histo_cov):
    cov = [0]
    cov_area2 = []
    for covered_position, count in cov_area.items():
        if count > 1:
            cov_area2.append(covered_position)

    cov_area = sorted(cov_area2, reverse=reverse)
    counter = 0
    for base_f in cov_area:
        count = False
        if reverse and base_f < peak_center:
            count = True
        elif not reverse and base_f > peak_center:
            count = True
        if count:
            if counter <= 3:
                counter += 1
                base_f = myround(base_f)
                if base_f in histo_cov[chrom]:
                    cov.append(histo_cov[chrom][base_f])
            else:
                break
    cov = max(cov)
    return cov, cov_area


def myround(x, base=10):
    '''Rounds to the nearest base'''
    return int(base * round(float(x) / base))


def find_peaks(density_dict, reverse, cutoff, histo_cov, side, peak_areas, chrom, CigarDict,start,end,splice_site_width,minimum_read_count,junctions):
    dist_range=[0]
    toWrite=[]
    for shift in range(1,splice_site_width+1):
        dist_range.append(shift)
        dist_range.append(-shift)

    entry_list = []
    for entry,density in density_dict.items():
        if len(density)>=minimum_read_count:
            entry_list.append([entry, density])
    for entry, density in sorted(entry_list, key=lambda x: len(x[1]), reverse=True):
        if entry in peak_areas[chrom][side]:
            continue
        best_extra_list, peak_center, cov_area, best_dirs,best_names = scan_for_best_bin(entry, dist_range, density_dict, peak_areas, chrom, side)
        cov, cov_area = determine_cov(cov_area, chrom, reverse, peak_center, histo_cov)

        if cov > 0:
            proportion = round(best_extra_list / cov, 3)
            if proportion > cutoff:
                Left_to_Right = best_dirs['+']
                Right_to_Left = best_dirs['-']
                Type = ''

                if Left_to_Right < Right_to_Left:
                    Type = '3' if reverse else '5'

                elif Left_to_Right > Right_to_Left:
                    Type = '5' if reverse else '3'

                if not Type:
                    continue
                passed = characterize_splicing_event([chrom, peak_center - splice_site_width, peak_center + splice_site_width], best_names, CigarDict,junctions)
                if passed:

                    start=str(peak_center - splice_site_width)
                    end=str(peak_center + splice_site_width)
                    toWrite.append([chrom,start,end,Type,side,str(proportion)])

                    for base in range(int(start), int(end) + 1):
                        peak_areas[chrom][side][base] = 1
        else:
            continue
    return peak_areas, toWrite


def collect_reads(reads, target_chrom,readAccuracy):

    histo_cov, histo_left_bases, histo_right_bases = {}, {}, {}
    histo_cov[target_chrom] = {}
    histo_left_bases[target_chrom] = {}
    histo_right_bases[target_chrom] = {}

    for line in open(reads):
        a = line.strip().split('\t')
        chrom = a[13]
        if chrom==target_chrom:
            name = a[9]
            dirn = a[8]
            length = int(a[10])
            begin, span = int(a[15]), int(a[16])
            blocksizes = a[18].split(',')[:-1]
            blockstarts = a[20].split(',')[:-1]
            cov_set = set()
            low_bounds, up_bounds = [], []
            aligned_bases = 0
            accuracy=readAccuracy[name]
            for x in range(0, len(blocksizes)):
                blockstart = int(blockstarts[x])
                blocksize = int(blocksizes[x])
                aligned_bases += blocksize
                blockend = blockstart + blocksize
                for y in range(0, blocksize, 10):
                    rounded = myround(blockstart + y)
                    cov_set.add(rounded)
                for yy in range(y, blocksize):
                    rounded = myround(blockstart + yy)
                    cov_set.add(rounded)
                if blockstart != begin:
                    up_bounds.append(blockstart)
                if blockend != span:
                    low_bounds.append(blockend)

            for rounded in cov_set:
                if rounded not in histo_cov[chrom]:
                    histo_cov[chrom][rounded] = 0
                histo_cov[chrom][rounded] += 1

            if accuracy<0.9:
                continue
            for low_bound in low_bounds:
                if low_bound not in histo_left_bases[chrom]:
                    histo_left_bases[chrom][low_bound] = []
                histo_left_bases[chrom][low_bound].append([name, begin, span, cov_set, dirn,accuracy,round(aligned_bases,-2)])
            for up_bound in up_bounds:
                if up_bound not in histo_right_bases[chrom]:
                    histo_right_bases[chrom][up_bound] = []
                histo_right_bases[chrom][up_bound].append([name, begin, span, cov_set, dirn,accuracy,round(aligned_bases,-2)])
    return histo_left_bases, histo_right_bases, histo_cov


def parse_genome(input_file, left_bounds, right_bounds, white_list_polyA):
    polyAWhiteList=[]
    chrom_list = set()
    gene_dict = {}
    for line in open(input_file):
        pAwl=False
        for element in white_list_polyA:
            if element in line:
                pAwl=True
        a = line.strip().split('\t')
        if len(a) <= 7:
            continue
        if a[2] == 'exon':
            testKey = a[8].split('transcript_id "')[1].split('"')[0]
            if not gene_dict.get(testKey):
                gene_dict[testKey] = []
            gene_dict[testKey].append((a[0], a[3], a[4], a[6], pAwl))

    for transcript_id in gene_dict:
        transcript_data = gene_dict[transcript_id]

        chrom = transcript_data[0][0]
        direction=transcript_data[0][3]
        chrom_list.add(chrom)
        if chrom not in right_bounds:
            left_bounds[chrom] = {'5': [], '3': []}
            right_bounds[chrom] = {'5': [], '3': []}

        start = sorted(transcript_data, key=lambda x: int(x[1]))[0][1]
        end = sorted(transcript_data, key=lambda x: int(x[2]), reverse=True)[0][2]
        direction=transcript_data[0][3]
        pAwl=transcript_data[0][4]
        if pAwl:
            if direction=='+':
                polyAWhiteList.append((chrom,direction,end,transcript_id))
            elif direction=='-':
                polyAWhiteList.append((chrom,direction,start,transcript_id))

        for entry in transcript_data:
            if entry[1] != start:
                if entry[3] == '+':
                    right_bounds[chrom]['3'].append(int(entry[1]) - 1)
                elif entry[3] == '-':
                    right_bounds[chrom]['5'].append(int(entry[1]) - 1)
            if entry[2] != end:
                if entry[3] == '+':
                    left_bounds[chrom]['5'].append(int(entry[2]))
                if entry[3] == '-':
                    left_bounds[chrom]['3'].append(int(entry[2]))
    return chrom_list, left_bounds, right_bounds, polyAWhiteList


def make_genome_bins(bounds, side, chrom, peak_areas,splice_site_width):
    toWrite=[]
    for type1 in ['5', '3']:
        covered = {}
        position_list = sorted(bounds[type1], key=int)
        for index1 in range(0, len(position_list)):
            if index1 in covered:
                continue
            sub_list = [position_list[index1]]
            for index2 in range(index1, len(position_list)):
                if position_list[index2] - max(sub_list) <= splice_site_width:
                    sub_list.append(position_list[index2])
                    covered[index2] = 1
                else:
                    break
            single = False
            if len(sub_list) > 1:
                splice_dists = []
                for splice_pos in range(0, len(sub_list) - 1):
                    splice_dists.append(int(sub_list[splice_pos + 1]) - int(sub_list[splice_pos]))
                if min(splice_dists) > 3:
                    for x in range(0, len(sub_list), 1):
                        if x != 0:
                            start = int(sub_list[x] - ((sub_list[x] - sub_list[x - 1]) / 2))
                        else:
                            start = int(sub_list[x]) - splice_site_width
                        if x != len(sub_list) - 1:
                            end = int(sub_list[x] + ((sub_list[x + 1] - sub_list[x]) / 2))
                        else:
                            end = int(sub_list[x]) + splice_site_width

                        toWrite.append([chrom,str(start),str(end),type1,side,'A'])
                        for base in range(start, end + 1):
                            peak_areas[chrom][side][base] = 1

                else:
                    single = True
            else:
                single = True
            if single:
                start = min(sub_list) - splice_site_width
                end = max(sub_list) + splice_site_width
                toWrite.append([chrom,str(start),str(end),type1,side,'A'])
                for base in range(start, end + 1):
                    peak_areas[chrom][side][base] = 1

    return peak_areas,toWrite



def get_chromosomes(infile,samfile,chrom_list,out_tmp,fastaList):
    outDict={}
    reads=[]
    new=False
    previous_chrom=''
    previous_start=0
    previous_end=0
    roots=set()
    print('\t\tparsing psl and splitting into loci')
    for line in open(infile):
        a=line.strip().split('\t')
        chrom=a[13]
        chrom_list.add(chrom)
        start=int(a[15])
        end=int(a[16])
        if chrom!=previous_chrom:
            new=True
        else:
            if start>previous_end:
                new=True
            else:
                previous_end=max(end,previous_end)
                reads.append(line)
        if new:
            if reads:
                root=previous_chrom+'~'+str(previous_start)+'~'+str(previous_end)
                fh=open(out_tmp+'/'+root+'.sam', 'w')
                fh.close()
                fh=open(out_tmp+'/'+root+'.fasta', 'w')
                fh.close()
                fh=open(out_tmp+'/'+root+'.psl', 'w')
                roots.add(root)

                for line in reads:
                    fh.write(line)
                    name=line.strip().split('\t')[9]
                    outDict[name]=root
                fh.close()
                reads=[]
            previous_chrom=chrom
            previous_end=end
            previous_start=start
            new=False
    if reads:
        root=previous_chrom+'~'+str(previous_start)+'~'+str(previous_end)
        fh=open(out_tmp+'/'+root+'.sam', 'w')
        fh.close()
        fh=open(out_tmp+'/'+root+'.fasta', 'w')
        fh.close()
        fh=open(out_tmp+'/'+root+'.psl', 'w')
        roots.add(root)
        for line in reads:
            fh.write(line)
            name=line.strip().split('\t')[9]
            outDict[name]=root
        fh.close()



    print('\t\tsplitting sam file into loci')
    used=set()
    for line in open(samfile):
        if line[0]=='@':
            continue
        a=line.strip().split('\t')
        name=a[0]
        if name not in used:
            if name in outDict:
                root=outDict[name]
                fh=open(out_tmp+'/'+root+'.sam', 'a')
                fh.write(line)
                fh.close()
            used.add(name)

    print('\t\tsplitting fasta file into loci')
    for fastaFile in fastaList:
        for name,seq,q in mp.fastx_read(fastaFile):
            if name in outDict:
                root=outDict[name]
                fh=open(out_tmp+'/'+root+'.fasta', 'a')
                fh.write('>%s\n%s\n' % (name,seq))
                fh.close()



def characterize_splicing_event(a,names,CigarDict,junctions):
    passed = False
    chromosome = a[0]
    splice_left = int(a[1])
    splice_right = int(a[2])
    readDict = []
    CIGARs=[]
    indeces = np.random.choice(np.arange(0, len(names)),
                                   min(len(names), 500), replace=False)
    for index in indeces:
        name = names[index]
        CIGARs.append(CigarDict[name])
    basecontext={}
    basecontext['allowed']=0
    basecontext['all']=0
    leftCS={}
    leftCS['Total']=0
    rightCS={}
    rightCS['Total']=0
    for status in ['*','+','-','=','|']:
        leftCS[status]=0
        rightCS[status]=0
    for CIGAR in CIGARs:
        bases,left,right = getCSaroundSS(CIGAR,splice_left,splice_right)
        basecontext['all'] += 1
#        if bases == 'gtag' or bases == 'ctac':
        if bases in junctions:
            basecontext['allowed'] += 1
        else:
            if bases not in basecontext:
                basecontext[bases]=0
            basecontext[bases]+=1
        if left:
            for thing in left:
                leftstatus = thing[0]
                leftCS[leftstatus]+=1
                leftCS['Total']+=1
        if right:
            for thing in right:
                rightstatus = thing[0]
                rightCS[rightstatus]+=1
                rightCS['Total']+=1

    basepass = False
    if basecontext['allowed']/basecontext['all']>0.85:
        basepass = True
        leftacc=leftCS['=']/leftCS['Total']
        rightacc=rightCS['=']/rightCS['Total']
    if basepass:
        if leftacc>0.85 and rightacc>0.85:
            passed = True
    CIGARs, CigarDict = 0, 0
    return passed



def find_ends(starts, ends, identity, count_dict, upstream_buffer,downstream_buffer,minimum_feature_count):
    start_peaks, end_peaks, extend = {}, {}, True
    if '+' in identity:
        count_dict['+'].add(identity)
    else:
        count_dict['-'].add(identity)


    start_count = {}
    end_count = {}
    for position in starts:
        if position not in start_count:
            start_count[position] = 0
        start_count[position] += 1

    for position in ends:
        if position not in end_count:
            end_count[position] = 0
        end_count[position] += 1

    show=False

    for position in sorted(starts):
        if position - upstream_buffer not in start_peaks:
            window_count = 0
            for i in np.arange(0,10):
                window_position = position + i
                if window_position in start_count:
                    window_count += start_count[window_position]
            if window_count >= minimum_feature_count:
                original_bin=[]
                for shift in np.arange(-upstream_buffer, downstream_buffer):
                    start_peaks[position + shift] = position
                    original_bin.append(position + shift)
                if extend:
                    bins=[]
                    for i in np.arange(min(original_bin),max(original_bin)):
                        bin=0
                        for shift in np.arange(0,10):
                            if i+shift in start_count:
                                 bin+=start_count[i+shift]
                        bins.append(bin)

                    best_bin=max(bins)
                    extended=True
                    adjacent=position - upstream_buffer
                    while extended:
                        window_count=0
                        adjacent_list=[]
                        for i in np.arange(1,11):
                            adjacent_pos = adjacent - i
                            adjacent_list.append(adjacent_pos)
                            if adjacent_pos in start_count:
                                window_count += start_count[adjacent_pos]
                        if best_bin > window_count >= minimum_feature_count:
                            for element in adjacent_list:
                                if element not in start_peaks:
                                    start_peaks[element] = position
                                    original_bin.append(element)
                                else:
                                    extended=False
                        else:
                            extended=False
                        adjacent = adjacent_pos
                        if extended:
                            count_dict['start_left'].add(identity)

                    extended = True
                    adjacent = position + downstream_buffer - 1
                    while extended:
                        window_count=0
                        adjacent_list=[]
                        for i in np.arange(1,11):
                            adjacent_pos = adjacent + i
                            adjacent_list.append(adjacent_pos)
                            if adjacent_pos in start_count:
                                window_count += start_count[adjacent_pos]
                        if best_bin > window_count >= minimum_feature_count:
                            for element in adjacent_list:
                                if element not in start_peaks:
                                    start_peaks[element] = position
                                    original_bin.append(element)
                                else:
                                     extended=False
                        else:
                            extended=False
                        adjacent = adjacent_pos
                        if extended:
                            count_dict['start_right'].add(identity)
    for position in sorted(ends, reverse=True):
        if position + upstream_buffer -1 not in end_peaks:
            window_count = 0
            for i in np.arange(0,10):
                window_position = position - i
                if window_position in end_count:
                    window_count += end_count[window_position]
            if window_count >= minimum_feature_count:
                original_bin=[]
                for shift in np.arange(-downstream_buffer, upstream_buffer):
                    end_peaks[position + shift] = position
                    original_bin.append(position + shift)
                if extend:
                    bins=[]
                    for i in np.arange(min(original_bin),max(original_bin)):
                        bin=0
                        for shift in np.arange(0,10):
                            if i+shift in end_count:
                                 bin+=end_count[i+shift]
                        bins.append(bin)

                    best_bin=max(bins)
                    extended = True
                    adjacent = position - downstream_buffer
                    while extended:
                        window_count=0
                        adjacent_list=[]
                        for i in np.arange(1,11):
                            adjacent_pos = adjacent - i
                            adjacent_list.append(adjacent_pos)
                            if adjacent_pos in end_count:
                                window_count += end_count[adjacent_pos]
                        if best_bin > window_count >= minimum_feature_count:
                            for element in adjacent_list:
                                if element not in end_peaks:
                                    end_peaks[element] = position
                                    original_bin.append(element)
                                else:
                                    extended=False
                        else:
                            extended=False
                        adjacent=adjacent_pos
                        if extended:
                            count_dict['end_left'].add(identity)

                    extended=True
                    adjacent=position + upstream_buffer - 1
                    while extended:
                        window_count=0
                        adjacent_list=[]
                        for i in np.arange(1,11):
                            adjacent_pos = adjacent + i
                            adjacent_list.append(adjacent_pos)
                            if adjacent_pos in end_count:
                                window_count += end_count[adjacent_pos]
                        if best_bin > window_count >= minimum_feature_count:
                            for element in adjacent_list:
                                if element not in end_peaks:
                                    original_bin.append(element)
                                    end_peaks[element] = position
                                else:
                                    extended=False
                        else:
                            extended=False
                        adjacent = adjacent_pos
                        if extended:
                            count_dict['end_right'].add(identity)

    return start_peaks, end_peaks, count_dict


def sort_reads_into_splice_junctions(splice_dict, fasta_file, infile):
    start_end_dict, start_end_dict_mono, readDict = {}, {}, {}
    for name,seq,q in mp.fastx_read(fasta_file):
        readDict[name] = [name, seq]

    for line in open(infile):
        a = line.strip().split('\t')
        read_chromosome, read_direction = a[13], a[8]
        name = a[9]
        show = False
        read_direction = '+'  # ignores read direction
        start, end = int(a[15]), int(a[16])
        if read_direction == '+':
            left_extra, right_extra = int(a[11]), int(a[10]) - int(a[12])
        if read_direction == '-':
            right_extra, left_extra = int(a[11]), int(a[10]) - int(a[12])

        failed = False
        identity = read_chromosome + '_'

        blocksizes = a[18].split(',')[:-1]
        blockstarts = a[20].split(',')[:-1]

        for x in range(0, len(blocksizes) - 1):
            blockstart = int(blockstarts[x])
            blocksize = int(blocksizes[x])
            left_splice = blockstart + blocksize
            right_splice = int(blockstarts[x + 1])
            if right_splice - left_splice > 50:
                if read_chromosome not in splice_dict:
                    failed = True
                    break
                else:
                    if not splice_dict[read_chromosome].get(left_splice) or not splice_dict[read_chromosome].get(right_splice):
                        failed = True
                        break
                left_splice_site = splice_dict[read_chromosome][left_splice]
                right_splice_site = splice_dict[read_chromosome][right_splice]
                identity += str(left_splice_site) + '-' + str(right_splice_site) + '~'
        if not failed:
            if identity.split('_')[1] != '':
                if identity not in start_end_dict:
                    start_end_dict[identity] = []
                start_end_dict[identity].append((start, end,
                                                 readDict[name],
                                                 left_extra,
                                                 right_extra,
                                                 read_direction))
            else:
                if identity not in start_end_dict_mono:
                    start_end_dict_mono[identity] = []

                start_end_dict_mono[identity].append((start, end,
                                                      readDict[name],
                                                      left_extra,
                                                      right_extra,
                                                      read_direction))
    return start_end_dict, start_end_dict_mono


def group_mono_exon_transcripts(start_end_dict, start_end_dict_mono):
    for identity in start_end_dict_mono:
        previous_end = 0
        iso_counter = 0
        new_identity = identity + 'M' + str(iso_counter)
        positions = start_end_dict_mono[identity]
        for position in sorted(positions):
            start = position[0]
            end = position[1]
            if start > previous_end:
                iso_counter += 1
                new_identity = identity + 'M' + str(iso_counter)
                if new_identity not in start_end_dict:
                    start_end_dict[new_identity] = []
                start_end_dict[new_identity].append(position)
                previous_end = max(end, previous_end)
            else:
                if new_identity not in start_end_dict:
                    start_end_dict[new_identity] = []
                start_end_dict[new_identity].append(position)
                previous_end = end

    return start_end_dict


def define_start_end_sites(start_end_dict, start_end_dict_mono,upstream_buffer,downstream_buffer,minimum_feature_count):
    left_extras, right_extras = {}, {}
    isoform_counter, isoform_dict, IsoDict = 0, {}, {}

    start_end_dict = group_mono_exon_transcripts(start_end_dict, start_end_dict_mono)
    number_of_isoforms = len(start_end_dict)
    counter = 0
    count_dict={}
    count_dict['start_left']=set()
    count_dict['start_right']=set()
    count_dict['end_left']=set()
    count_dict['end_right']=set()
    count_dict['+']=set()
    count_dict['-']=set()

    for identity in sorted(start_end_dict):
        counter += 1
        starts, ends = [], []
        ends = []
        positions = start_end_dict[identity]
        length_of_positions = len(positions)
        indexes = np.random.choice(range(0, length_of_positions, 1),
                                   size=min(10000, length_of_positions),
                                   replace=False)
        subsampled_positions = []
        for index in indexes:
            subsampled_positions.append(positions[index])
        for position in subsampled_positions:
            starts.append(int(position[0]))
            ends.append(int(position[1]))


        start_dict, end_dict, count_dict = find_ends(starts, ends, identity, count_dict,upstream_buffer,downstream_buffer,minimum_feature_count)
        matched_positions = []
        left_extras[identity], right_extras[identity] = {}, {}
        for start, end, read, left_extra, right_extra, read_direction in positions:
            if int(start) in start_dict and int(end) in end_dict:
                left = start_dict[int(start)]
                right = end_dict[int(end)]
                if not left_extras[identity].get((left, right)):
                    left_extras[identity][(left, right)] = []
                    right_extras[identity][(left, right)] = []

                left_extras[identity][(left, right)].append(int(left_extra))
                right_extras[identity][(left, right)].append(int(right_extra))
                matched_positions.append((left, right, read, read_direction))

        median_left = {}
        median_right = {}
        for combination, values in left_extras[identity].items():
            median = np.median(values)
            median_left[combination] = median
        for combination, values in right_extras[identity].items():
            median = np.median(values)
            median_right[combination] = median

        for left, right, read, read_direction in matched_positions:
            medianLeft = median_left[(left, right)]
            medianRight = median_right[(left, right)]
            new_identity = identity + '_' + read_direction + '_' + str(left) \
                + '_' + str(right) + '_' \
                + str(round(medianLeft, 2)) \
                + '_' + str(round(medianRight, 2))
            if new_identity not in isoform_dict:
                isoform_counter += 1
                isoform_dict[new_identity] = isoform_counter

            IsoformName = str(isoform_dict[new_identity])
            if IsoformName not in IsoDict:
                IsoDict[IsoformName]=[]
            IsoDict[IsoformName].append(read)
    return IsoDict

def revComp(seq):
    '''Returns the reverse complement of a seq'''
    bases = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', '-': '-'}
    return ''.join([bases[x] for x in list(seq)])[::-1]


def determine_consensus(reads):
    '''Aligns and returns the consensus'''
    corrected_consensus = ''
    repeats = '0'
    fasta_reads = []
    names=[]
    for read, seq in reads:
        fasta_reads.append((read, seq))
        names.append(read)
    repeats = str(len(fasta_reads))

    indeces = np.random.choice(np.arange(0, len(fasta_reads)),
                                   min(len(fasta_reads), 100), replace=False)
    subsample_fasta_reads = []
    for index in indeces:
        subsample_fasta_reads.append(fasta_reads[index])
    first = subsample_fasta_reads[0][1]
    sequences=[]
    seq_lengths=[]
    mm_align = mp.Aligner(seq=first, preset='map-ont')
    for read,sequence in subsample_fasta_reads:
        seq_lengths.append(len(sequence))
        for hit in mm_align.map(sequence):
             if hit.is_primary:
                 if hit.strand==1:
                     sequences.append(sequence)
                 elif hit.strand==-1:
                     sequences.append(mp.revcomp(sequence))

    if np.median(seq_lengths)>20000:
#        consensus_sequence=split_consensus(sequences)
        consensus_sequence = sequences[0]
    else:
        res = poa_aligner.msa(sequences, out_cons=True, out_msa=False)
        if len(sequences)<=2:
            consensus_sequence = sequences[0]
        elif not res.cons_seq:
            consensus_sequence = sequences[0]
        else:
            consensus_sequence = res.cons_seq[0]
    res, mm_align, sequences = 0, 0, 0
    return consensus_sequence,names


#def split_consensus(sequences):

